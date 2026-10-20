# Regression tests for inference changes. Run standalone after `using Test,
# Random, LinearAlgebra, Pathogen`, or through runtests.jl.
perf_one(p, pop, i) = 1.0
perf_constant(p, pop, i) = p[1]
perf_pair(p, pop, i, j) = p[1] * (1.0 + abs(i-j) / 10)

@testset "Sparse sampled networks" begin
  n = 1000
  small, large = TransmissionNetwork(n), TransmissionNetwork(2n)
  @test size(small.internal) == (n, n)
  @test Base.summarysize(large) < 3 * Base.summarysize(small)
  Pathogen.update!(small, ExogenousTransmission(1))
  for i in 2:n
    Pathogen.update!(small, EndogenousTransmission(i, i-1))
  end
  @test sum(small) == n
  @test Pathogen.nnz(small.internal) == n-1
  copied = copy(small)
  copied.internal[1, 2] = false
  @test small.internal[1, 2]
  @test sum(innerjoin(small, small)) == n
  @test sum(leftantijoin(small, small)) == 0
  dense = TransmissionNetwork(copy(small.external), BitMatrix(small.internal))
  @test dense.internal == small.internal
  tnd = TNDistribution([small, dense])
  @test sum(tnd) == n
  @test mode(tnd).internal == small.internal
  bad = falses(3, 3)
  bad[1, 3] = bad[2, 3] = true
  @test_throws ErrorException TransmissionNetwork(falses(3), bad)
end

function perf_fixture(::Type{S}; n=8, phylo=false) where S
  pop = Population(DataFrame(x=collect(1:n)))
  funcs = Any[perf_constant, perf_one, perf_pair, perf_one]
  params = [[0.05], Float64[], [0.3], Float64[]]
  priors = [UnivariateDistribution[Uniform(0.01, 0.1)], UnivariateDistribution[],
            UnivariateDistribution[Uniform(0.1, 0.8)], UnivariateDistribution[]]
  for state in (State_E, State_R)
    if state ∈ S
      push!(funcs, perf_constant)
      push!(params, [0.5])
      push!(priors, UnivariateDistribution[Uniform(0.1, 1.0)])
    end
  end
  rf, rp = RiskFunctions{S}(funcs...), RiskParameters{S}(params...)
  states = [State_I; fill(State_S, n-1)]
  sim = phylo ? Simulation(pop, states, rf, rp, JC69(0.05)) : Simulation(pop, states, rf, rp)
  simulate!(MersenneTwister(42), sim; tmax=100.0)
  obs = State_R ∈ S ? observe(MersenneTwister(43), sim, Uniform(0.0, 0.1), Uniform(0.0, 0.1); seq_len=8, force=true) :
                     observe(MersenneTwister(43), sim, Uniform(0.0, 0.1); seq_len=8)
  ee = S == SEIR ? EventExtents{S}(30.0, 0.2, 0.2) : S == SEI ? EventExtents{S}(30.0, 0.2) :
       S == SIR ? EventExtents{S}(0.2, 0.2) : EventExtents{S}(0.2)
  mcmc = phylo ? MCMC(obs, ee, pop, states, rf, RiskPriors{S}(priors...), JC69abs, UnivariateDistribution[Uniform(0.01, 0.1)]) :
                 MCMC(obs, ee, pop, states, rf, RiskPriors{S}(priors...))
  start!(MersenneTwister(44), mcmc; attempts=100)
  return sim, mcmc
end

@testset "Incremental removals preserve hazards and historical weights" begin
  for S in (SIR, SEIR)
    sim, _ = perf_fixture(S)
    states = copy([State_I; fill(State_S, Pathogen.individuals(sim.population)-1)])
    tr = Pathogen.initialize(Pathogen.TransmissionRates, states, sim.population, sim.risk_functions, sim.risk_parameters)
    rates = Pathogen.initialize(Pathogen.EventRates, tr, states, sim.population, sim.risk_functions, sim.risk_parameters)
    history = [(0.1, 2, State_I), (0.2, 1, State_R), (0.3, 3, State_I), (0.4, 2, State_R), (0.5, 3, State_R)]
    for (t, id, state) in history
      if S == SEIR && state == State_I
        exposure = Pathogen.Event{S}(t-0.01, id, State_E)
        Pathogen.update!(states, exposure)
        Pathogen.update!(tr, rates, exposure, states, sim.population, sim.risk_functions, sim.risk_parameters)
      end
      historical = copy(tr.internal[:, states .!= Ref(State_S)])
      event = Pathogen.Event{S}(t, id, state)
      Pathogen.update!(states, event)
      Pathogen.update!(tr, rates, event, states, sim.population, sim.risk_functions, sim.risk_parameters)
      expected_tr = Pathogen.initialize(Pathogen.TransmissionRates, states, sim.population, sim.risk_functions, sim.risk_parameters)
      expected = Pathogen.initialize(Pathogen.EventRates, expected_tr, states, sim.population, sim.risk_functions, sim.risk_parameters)
      @test rates.infection ≈ expected.infection
      @test rates.removal ≈ expected.removal
      if S == SEIR
        @test rates.exposure ≈ expected.exposure
      end
      if state == State_R
        @test tr.internal[:, states .!= Ref(State_S)] == historical
      end
    end
  end
end

@testset "Complete likelihood before rejection" begin
  pop = Population(DataFrame(x=[1]))
  rf = RiskFunctions{SIR}(perf_constant, perf_one, perf_pair, perf_one, perf_constant)
  rp = RiskParameters{SIR}([1.0], Float64[], [0.0], Float64[], [1000.0])
  args = (rp, rf, Events{SIR}([0.0], [0.001]), pop, [State_S])
  full = Pathogen.loglikelihood(args...)[1]
  @test full ≈ log(1000.0) - 1.0
  @test Pathogen.loglikelihood(args...; early_decision_value=1.0)[1] == full
end

@testset "Cache follows accepted parameters across iterations" begin
  for phylo in (false, true)
    _, mcmc = perf_fixture(SIR; phylo=phylo)
    cached, fresh = deepcopy(mcmc.markov_chains[1]), deepcopy(mcmc.markov_chains[1])
    a, b = MersenneTwister(45), MersenneTwister(45)
    cache = Pathogen.TransmissionRateCache(Pathogen.individuals(mcmc.population))
    Σ = Matrix{Float64}(I, length(cached.risk_parameters[end]), length(cached.risk_parameters[end])) * 0.01
    for _ in 1:16
      if phylo
        Pathogen.update!(a, cached, mcmc, Σ, fill(0.0001, 1, 1), 0.02; transmission_rate_cache=cache)
        Pathogen.update!(b, fresh, mcmc, Σ, fill(0.0001, 1, 1), 0.02)
      else
        Pathogen.update!(a, cached, mcmc, Σ, 0.02; transmission_rate_cache=cache)
        Pathogen.update!(b, fresh, mcmc, Σ, 0.02)
      end
      @test cached.log_posterior[end] ≈ fresh.log_posterior[end]
      @test convert(Vector{Float64}, cached.risk_parameters[end]) == convert(Vector{Float64}, fresh.risk_parameters[end])
      @test isequal(convert(Vector{Float64}, cached.events[end]), convert(Vector{Float64}, fresh.events[end]))
      @test cached.transmission_network[end].internal == fresh.transmission_network[end].internal
      @test rand(copy(a)) == rand(copy(b))
    end
  end
end
