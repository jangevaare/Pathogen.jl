# Regression tests for inference changes. Run standalone after `using Test,
# Random, LinearAlgebra, Pathogen`, or through runtests.jl.
perf_one(p, pop, i) = 1.0
perf_constant(p, pop, i) = p[1]
perf_pair(p, pop, i, j) = p[1] * (1.0 + abs(i-j) / 10)
perf_distance_pair(p, pop, i, j) = pop.distances[j, i]
perf_covariate_rate(p, pop, i) = pop.risks[i, :rate]

@testset "Dominant-rate removal preserves small surviving hazards" begin
  rf = RiskFunctions{SIR}(perf_constant, perf_one, perf_distance_pair, perf_one, perf_covariate_rate)
  rp = RiskParameters{SIR}([0.0], Float64[], Float64[], Float64[], Float64[])
  pop = Population(DataFrame(rate=ones(3)), [0.0 0.0 1e16; 0.0 0.0 1.0; 0.0 0.0 0.0])
  states = [State_I, State_I, State_S]
  tr = Pathogen.initialize(Pathogen.TransmissionRates, states, pop, rf, rp)
  rates = Pathogen.initialize(Pathogen.EventRates, tr, states, pop, rf, rp)
  total = Ref(sum(rates))
  event = Pathogen.Event{SIR}(0.0, 1, State_R)
  Pathogen.update!(states, event)
  Pathogen.update!(tr, rates, event, states, pop, rf, rp; rate_total=total)
  @test rates.infection[3] == 1.0
  @test total[] == sum(rates) == 2.0
  events = Events{SIR}([-Inf, -Inf, 0.01], [0.0, NaN, NaN])
  ll = Pathogen.loglikelihood(rp, rf, events, pop, [State_I, State_I, State_S])[1]
  @test ll ≈ log(1.0 / (1e16 + 2.0)) - 0.02

  pop = Population(DataFrame(rate=[1e16, 1.0, 2.0]), zeros(3, 3))
  events = Events{SIR}(fill(-Inf, 3), [0.0, 1.0, NaN])
  work = Pathogen.LikelihoodWorkspace(SIR, 3)
  ll = Pathogen.loglikelihood(rp, rf, events, pop, fill(State_I, 3); workspace=work)[1]
  @test ll ≈ -3.0
  @test work.rate_total[] == sum(work.event_rates) == 2.0
end

@testset "Child traversal and significant transmission propagation" begin
  for shape in (:chain, :star, :forest)
    n = 12
    network = TransmissionNetwork(n)
    network.external[1] = true
    for i in 2:n
      if shape == :forest && i == 7
        network.external[i] = true
      else
        parent = shape == :chain ? i-1 : shape == :star ? 1 : (i < 7 ? i-1 : 7)
        Pathogen.update!(network, EndogenousTransmission(i, parent))
      end
    end
    children = Pathogen._children(network)
    @test sum(length, children) == n - sum(network.external)
    @test Pathogen._pathway_from(1, network; depth=1) == [1; children[1]]
    @test Pathogen._source(1, network) === nothing
    for i in 2:n
      @test Pathogen._source(i, network) == findfirst(network.internal[:, i])
    end
    times = [collect(0.0:n-1) fill(NaN, n)]
    times[network.external, 1] .= -Inf
    for observed in (Int[], [n], [2, 5, n])
      times[:, 2] .= NaN
      times[observed, 2] = observed .+ 1.0
      expected = [times[i, 1] == -Inf || any(times[Pathogen._pathway_from(i, network), 2] .>= times[i, 1]) for i in 1:n]
      @test Pathogen._significant_transmissions(times, network, children) == expected
    end
  end
  cyclic = TransmissionNetwork(2)
  cyclic.internal[1, 2] = cyclic.internal[2, 1] = true
  @test_throws ArgumentError Pathogen._significant_transmissions(zeros(2, 2), cyclic, Pathogen._children(cyclic))
end

@testset "Event proposals own their vectors and cover every event" begin
  for S in (SI, SIR, SEI, SEIR)
    states = convert(DiseaseStates, S)[2:end]
    original = Events{S}((fill(Float64(i), 4) for i in eachindex(states))...)
    proposal = copy(original)
    Pathogen.update!(proposal, Pathogen.Event{S}(12.0, 2, states[1]))
    @test original[states[1]][2] == 1.0
    @test proposal[states[1]][2] == 12.0
    for state in states
      @test original[state] !== proposal[state]
    end
  end
  for n in (1, 2, 7, 13), batches in 1:n
    indices = reduce(vcat, [collect(Pathogen._event_batch(n, batches, b)) for b in 1:batches])
    @test indices == collect(1:n)
  end
  @test_throws DimensionMismatch Events{SEIR}(zeros(2), zeros(3), zeros(2))
  @test_throws DimensionMismatch Events{SEI}(zeros(2), zeros(3))
  @test_throws DimensionMismatch Events{SIR}(zeros(2), zeros(3))
end

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

# Deliberately slow oracle: initialize all current hazards from scratch before
# each event, independently of incremental updates and workspace reuse.
function perf_reference_ll(sim)
  S = typeof(sim.events).parameters[1]
  states = [State_I; fill(State_S, Pathogen.individuals(sim.population)-1)]
  times = convert(Vector{Float64}, sim.events)
  kinds = convert(DiseaseStates, S)[2:end]
  ll, last = 0.0, nothing
  n = length(states)
  for k in sortperm(times)
    t = times[k]
    t == -Inf && continue
    isnan(t) && break
    state, id = kinds[fld(k-1, n)+1], mod1(k, n)
    tr = Pathogen.initialize(Pathogen.TransmissionRates, states, sim.population, sim.risk_functions, sim.risk_parameters)
    rates = Pathogen.initialize(Pathogen.EventRates, tr, states, sim.population, sim.risk_functions, sim.risk_parameters)
    total = sum(rates)
    event_rate = rates[state][id]
    ll += last === nothing ? log(event_rate / total) : log(event_rate) - total * (t-last)
    states[id] = state
    last = t
  end
  return ll
end

@testset "Reusable likelihood matches independent hazard replay" begin
  for S in (SI, SIR, SEI, SEIR)
    sim, _ = perf_fixture(S)
    n = Pathogen.individuals(sim.population)
    args = (sim.risk_parameters, sim.risk_functions, sim.events, sim.population, [State_I; fill(State_S, n-1)])
    work = Pathogen.LikelihoodWorkspace(S, n)
    expected = perf_reference_ll(sim)
    ll, tr, tx = Pathogen.loglikelihood(args...; workspace=work)
    @test ll ≈ expected
    @test work.rate_total[] ≈ sum(work.event_rates) atol=1e-10
    @test Pathogen.loglikelihood(args...)[1] ≈ expected
    @test Pathogen.loglikelihood(args...; workspace=work)[1] ≈ expected
    # In-place edits of an existing Events object must invalidate cached ordering.
    for state in convert(DiseaseStates, S)[2:end], i in 1:n
      t = sim.events[state][i]
      isfinite(t) && (sim.events[state][i] += 0.1)
    end
    @test Pathogen.loglikelihood(args...; workspace=work)[1] ≈ perf_reference_ll(sim)
    @test work.event_order == sortperm(convert(Vector{Float64}, sim.events))
    @test Pathogen.loglikelihood(args...; workspace=work, loglikelihood_output=false,
                               transmission_rates_output=false, transmissions_output=false) == (nothing, nothing, nothing)
    @test tr === work.transmission_rates
    @test_throws DimensionMismatch Pathogen.loglikelihood(args...; workspace=Pathogen.LikelihoodWorkspace(S, n+1))
  end
  events = Events{SIR}([0.0, 1.0], [2.0, 3.0])
  work = Pathogen.LikelihoodWorkspace(SIR, 2)
  before = copy(Pathogen._event_order!(work, events))
  events.removal[1] = 0.5
  @test Pathogen._event_order!(work, events) != before
  @test work.event_order == sortperm(convert(Vector{Float64}, events))
end

@testset "Empty augmentation and explicit batch validation" begin
  _, mcmc = perf_fixture(SI; n=1)
  mc = mcmc.markov_chains[1]
  Σ = Matrix{Float64}(I, length(mc.risk_parameters[end]), length(mc.risk_parameters[end])) * 0.0001
  Pathogen.update!(MersenneTwister(1), mc, mcmc, Σ, 0.01; event_batches=0)
  @test mc.iterations == 1
  @test_throws ArgumentError Pathogen.update!(MersenneTwister(1), mc, mcmc, Σ, 0.01; event_batches=-1)
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
    work = Pathogen.LikelihoodWorkspace(SIR, Pathogen.individuals(mcmc.population))
    Σ = Matrix{Float64}(I, length(cached.risk_parameters[end]), length(cached.risk_parameters[end])) * 0.01
    for _ in 1:16
      if phylo
        Pathogen.update!(a, cached, mcmc, Σ, fill(0.0001, 1, 1), 0.02; transmission_rate_cache=cache, likelihood_workspace=work)
        Pathogen.update!(b, fresh, mcmc, Σ, fill(0.0001, 1, 1), 0.02)
      else
        Pathogen.update!(a, cached, mcmc, Σ, 0.02; transmission_rate_cache=cache, likelihood_workspace=work)
        Pathogen.update!(b, fresh, mcmc, Σ, 0.02)
      end
      @test work.proposal_cache.internal !== cache.internal
      @test cached.log_posterior[end] ≈ fresh.log_posterior[end]
      @test convert(Vector{Float64}, cached.risk_parameters[end]) == convert(Vector{Float64}, fresh.risk_parameters[end])
      @test isequal(convert(Vector{Float64}, cached.events[end]), convert(Vector{Float64}, fresh.events[end]))
      @test cached.transmission_network[end].internal == fresh.transmission_network[end].internal
      @test rand(copy(a)) == rand(copy(b))
    end
  end
end
