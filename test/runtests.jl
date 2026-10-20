using Test, Distributed, Random, LinearAlgebra, Distances, Pathogen, StableRNGs
import PhyloTrees._treeplot

addprocs(3)

@everywhere using Pathogen
@everywhere include(joinpath(@__DIR__, "risk_functions.jl"))

# Set RNG seed
rng = StableRNG(543210)

# using Logging
# global_logger(ConsoleLogger(stderr, LogLevel(-1000)))

# Define population
n = 25
risks = DataFrame(x = rand(Uniform(0, 5), n),
                  y = rand(Uniform(0, 5), n),
                  riskfactor1 = rand(Gamma(), n))

pop = Population(risks,
                 [euclidean([risks[i, :x]; risks[i, :y]], [risks[j, :x]; risks[j, :y]]) for i = 1:n, j = 1:n])

@testset "SEIR TNILM" begin
  # Some commonly used functions/examples provided in helpers/RiskFunctions.jl
  # For SEIR, risk functions and parameters in order of: sparks, susceptibility, infectivity, transmissibility, latency, and removal
  rf = RiskFunctions{SEIR}(_constant,
                           _coefficient,
                           _powerlaw,
                           _one,
                           _constant,
                           _constant)

  rparams = RiskParameters{SEIR}([0.01],
                                 [1.0],
                                 [2.0, 5.0],
                                 Float64[],
                                 [0.1],
                                 [0.05])

  sim = Simulation(pop,
                   rf,
                   rparams)

  simulate!(rng, sim, tmax=100.0)

  state_counts = [Pathogen._count_by_state(sim.events, State_S, 50.0)
                  Pathogen._count_by_state(sim.events, State_E, 50.0)
                  Pathogen._count_by_state(sim.events, State_I, 50.0)
                  Pathogen._count_by_state(sim.events, State_R, 50.0)]

  @test sum(state_counts) == n
  @test length(sim.disease_states) == n
  @test size(sim.transmission_network.internal) == (n, n)

  obs = observe(rng, sim, Uniform(), Uniform())
  @test length(obs.infection) == n

  rpriors = RiskPriors{SEIR}([Uniform(0.0, 0.1)],
                             [Uniform(0.0, 2.0)],
                             [Uniform(0.0, 4.0); Uniform(1.0, 8.0)],
                             UnivariateDistribution[],
                             [Uniform(0.0, 1.0)],
                             [Uniform(0.0, 1.0)])
  ee = EventExtents{SEIR}((2.0, 20.0), 2.0, 2.0)
  mcmc = MCMC(obs, ee, pop, rf, rpriors)
  start!(mcmc, 3, attempts=100)
  @test length(mcmc.markov_chains) == 3
  @test all([length(mcmc.markov_chains[i].risk_parameters[1]) for i=1:3] .== length(rparams))

  # Check bounds of Events initialization
  for j = 1:1, i = 1:n
    if mcmc.markov_chains[j].events[1].exposure[i] > -Inf
      l, u = Pathogen._bounds(i, State_E, ee, obs, mcmc.markov_chains[j].events[1], mcmc.markov_chains[j].transmission_network[1])
      @test l < u
    end
    if mcmc.markov_chains[j].events[1].infection[i] > -Inf
      l, u = Pathogen._bounds(i, State_I, ee, obs, mcmc.markov_chains[j].events[1], mcmc.markov_chains[j].transmission_network[1])
      @test l < u
    end
    if mcmc.markov_chains[j].events[1].removal[i] > -Inf
      l, u = Pathogen._bounds(i, State_R, ee, obs, mcmc.markov_chains[j].events[1], mcmc.markov_chains[j].transmission_network[1])
      @test l < u
    end
  end
  iterate!(mcmc, 100, 1.0, condition_on_network=true)
  @test all([mcmc.markov_chains[i].iterations for i = 1:3] .== 100)
  tnd = TNDistribution(mcmc, burnin=2)
  @test sum(tnd) ≈ sum(obs.infection .!== NaN)
  @test sum(mode(tnd)) ≈ sum(obs.infection .!== NaN)
  @test all(Pathogen._outdegree(tnd) .>= 0.0)
  mcmc = MCMC(obs, ee, pop, rf, rpriors, tnprior=tnd)
  start!(mcmc, attempts=100)
  iterate!(mcmc, 100, 1.0, condition_on_network=false)
  @test length(mcmc.markov_chains) == 1
  @test typeof(summary(mcmc, burnin=50, thin=2)) == DataFrame

  tn = mcmc.markov_chains[1].transmission_network[end]
  for i = 1:n
    pt = Pathogen._pathway_to(i, tn, depth=1)
    if pt != nothing
      @test i in Pathogen._pathway_from(pt[2], tn, depth=1)
    end
  end
end

@testset "SEI TNILM" begin
  rf = RiskFunctions{SEI}(_constant,
                          _coefficient,
                          _powerlaw,
                          _one,
                          _constant)

  rparams = RiskParameters{SEI}([0.01],
                                [1.0],
                                [2.0, 5.0],
                                Float64[],
                                [0.1])

  sim = Simulation(pop,
                   rf,
                   rparams)

  simulate!(rng, sim, tmax=100.0)

  state_counts = [Pathogen._count_by_state(sim.events, State_S, 50.0)
                  Pathogen._count_by_state(sim.events, State_E, 50.0)
                  Pathogen._count_by_state(sim.events, State_I, 50.0)]

  @test sum(state_counts) == n

  @test length(sim.disease_states) == n
  @test size(sim.transmission_network.internal) == (n, n)

  obs = observe(rng, sim, Uniform())
  @test length(obs.infection) == n

  rpriors = RiskPriors{SEI}([Uniform(0.0, 0.1)],
                            [Uniform(0.0, 2.0)],
                            [Uniform(0.0, 4.0); Uniform(1.0, 8.0)],
                            UnivariateDistribution[],
                            [Uniform(0.0, 1.0)])

  ee = EventExtents{SEI}((2.0, 20.0), 2.0)
  mcmc = MCMC(obs, ee, pop, rf, rpriors)
  start!(mcmc, 3, attempts = 100)
  @test length(mcmc.markov_chains) == 3
  @test all([length(mcmc.markov_chains[i].risk_parameters[1]) for i=1:3] .== length(rparams))

  # Check bounds of Events initialization
  for j = 1:3, i = 1:n
    if mcmc.markov_chains[j].events[1].exposure[i] > -Inf
      l, u = Pathogen._bounds(i, State_E, ee, obs, mcmc.markov_chains[j].events[1], mcmc.markov_chains[j].transmission_network[1])
      @test l < u
    end
    if mcmc.markov_chains[j].events[1].infection[i] > -Inf
      l, u = Pathogen._bounds(i, State_I, ee, obs, mcmc.markov_chains[j].events[1], mcmc.markov_chains[j].transmission_network[1])
      @test l < u
    end
  end
  iterate!(mcmc, 100, 1.0, condition_on_network=true)
  @test all([mcmc.markov_chains[i].iterations for i = 1:3] .== 100)
  tnd = TNDistribution(mcmc, burnin=2)
  @test sum(tnd) ≈ sum(obs.infection .!== NaN)
  @test sum(mode(tnd)) ≈ sum(obs.infection .!== NaN)
  @test all(Pathogen._outdegree(tnd) .>= 0.0)
  mcmc = MCMC(obs, ee, pop, rf, rpriors, tnprior=tnd)
  start!(mcmc, attempts=50)
  iterate!(mcmc, 100, 1.0, condition_on_network=false)
  @test length(mcmc.markov_chains) == 1
  @test typeof(summary(mcmc, burnin=50, thin=2)) == DataFrame

  tn = mcmc.markov_chains[1].transmission_network[end]
  for i = 1:n
    pt = Pathogen._pathway_to(i, tn, depth=1)
    if pt != nothing
      @test i in Pathogen._pathway_from(pt[2], tn, depth=1)
    end
  end
end


@testset "SIR TNILM" begin
  rf = RiskFunctions{SIR}(_constant,
                          _coefficient,
                          _powerlaw,
                          _one,
                          _constant)

  rparams = RiskParameters{SIR}([0.01],
                                [1.0],
                                [2.0, 5.0],
                                Float64[],
                                [0.05])

  sim = Simulation(pop,
                   rf,
                   rparams)

  simulate!(rng, sim, tmax=100.0)

  state_counts = [Pathogen._count_by_state(sim.events, State_S, 50.0)
                  Pathogen._count_by_state(sim.events, State_I, 50.0)
                  Pathogen._count_by_state(sim.events, State_R, 50.0)]

  @test sum(state_counts) == n

  @test length(sim.disease_states) == n
  @test size(sim.transmission_network.internal) == (n, n)

  obs = observe(rng, sim, Uniform(), Uniform())
  @test length(obs.infection) == n

  rpriors = RiskPriors{SIR}([Uniform(0.0, 0.1)],
                            [Uniform(0.0, 2.0)],
                            [Uniform(0.0, 4.0); Uniform(1.0, 8.0)],
                            UnivariateDistribution[],
                            [Uniform(0.0, 1.0)])

  ee = EventExtents{SIR}(2.0, 2.0)
  mcmc = MCMC(obs, ee, pop, rf, rpriors)
  start!(mcmc, 3, attempts=100)
  @test length(mcmc.markov_chains) == 3
  @test all([length(mcmc.markov_chains[i].risk_parameters[1]) for i=1:3] .== length(rparams))

  # Check bounds of Events initialization
  for j = 1:3, i = 1:n
    if mcmc.markov_chains[j].events[1].infection[i] > -Inf
      l, u = Pathogen._bounds(i, State_I, ee, obs, mcmc.markov_chains[j].events[1], mcmc.markov_chains[j].transmission_network[1])
      @test l < u
    end
    if mcmc.markov_chains[j].events[1].removal[i] > -Inf
      l, u = Pathogen._bounds(i, State_R, ee, obs, mcmc.markov_chains[j].events[1], mcmc.markov_chains[j].transmission_network[1])
      @test l < u
    end
  end
  iterate!(mcmc, 100, 1.0, condition_on_network=true)
  @test all([mcmc.markov_chains[i].iterations for i = 1:3] .== 100)
  tnd = TNDistribution(mcmc, burnin=2)
  @test sum(tnd) ≈ sum(obs.infection .!== NaN)
  @test sum(mode(tnd)) ≈ sum(obs.infection .!== NaN)
  @test all(Pathogen._outdegree(tnd) .>= 0.0)
  mcmc = MCMC(obs, ee, pop, rf, rpriors, tnprior=tnd)
  start!(mcmc, attempts=50)
  iterate!(mcmc, 100, 1.0, condition_on_network=false)
  @test length(mcmc.markov_chains) == 1
  @test typeof(summary(mcmc, burnin=50, thin=2)) == DataFrame

  tn = mcmc.markov_chains[1].transmission_network[end]
  for i = 1:n
    pt = Pathogen._pathway_to(i, tn, depth=1)
    if pt != nothing
      @test i in Pathogen._pathway_from(pt[2], tn, depth=1)
    end
  end
end

@testset "SI TNILM" begin
  rf = RiskFunctions{SI}(_constant,
                         _coefficient,
                         _powerlaw,
                         _one)

  rparams = RiskParameters{SI}([0.01],
                               [1.0],
                               [2.0, 5.0],
                               Float64[])

  sim = Simulation(pop,
                   rf,
                   rparams)

  simulate!(rng, sim, tmax=100.0)

  state_counts = [Pathogen._count_by_state(sim.events, State_S, 50.0)
                  Pathogen._count_by_state(sim.events, State_I, 50.0)]

  @test sum(state_counts) == n

  @test length(sim.disease_states) == n
  @test size(sim.transmission_network.internal) == (n, n)

  obs = observe(rng, sim, Uniform())
  @test length(obs.infection) == n

  rpriors = RiskPriors{SI}([Uniform(0.0, 0.1)],
                           [Uniform(0.0, 2.0)],
                           [Uniform(0.0, 4.0); Uniform(1.0, 8.0)],
                           UnivariateDistribution[])

  ee = EventExtents{SI}(2.0)
  mcmc = MCMC(obs, ee, pop, rf, rpriors)
  start!(mcmc, 3, attempts=100)
  @test length(mcmc.markov_chains) == 3
  @test all([length(mcmc.markov_chains[i].risk_parameters[1]) for i=1:3] .== length(rparams))
  # Check bounds of Events initialization
  for j = 1:3, i = 1:n
    if mcmc.markov_chains[j].events[1].infection[i] > -Inf
      l, u = Pathogen._bounds(i, State_I, ee, obs, mcmc.markov_chains[j].events[1], mcmc.markov_chains[j].transmission_network[1])
      @test l < u
    end
  end
  iterate!(mcmc, 100, 1.0, condition_on_network=true)
  @test all([mcmc.markov_chains[i].iterations for i = 1:3] .== 100)
  tnd = TNDistribution(mcmc, burnin=2)
  @test sum(tnd) ≈ sum(obs.infection .!== NaN)
  @test sum(mode(tnd)) ≈ sum(obs.infection .!== NaN)
  @test all(Pathogen._outdegree(tnd) .>= 0.0)
  mcmc = MCMC(obs, ee, pop, rf, rpriors, tnprior=tnd)
  start!(mcmc, attempts=50)
  iterate!(mcmc, 100, 1.0, condition_on_network=false)
  @test length(mcmc.markov_chains) == 1
  @test typeof(summary(mcmc, burnin=50, thin=2)) == DataFrame

  tn = mcmc.markov_chains[1].transmission_network[end]
  for i = 1:n
    pt = Pathogen._pathway_to(i, tn, depth=1)
    if pt != nothing
      @test i in Pathogen._pathway_from(pt[2], tn, depth=1)
    end
  end
end

@testset "SEIR PhyloILM" begin
  # Some commonly used functions/examples provided in helpers/RiskFunctions.jl
  # For SEIR, risk functions and parameters in order of: sparks, susceptibility, infectivity, transmissibility, latency, and removal
  rf = RiskFunctions{SEIR}(_zero,
                           _coefficient,
                           _powerlaw,
                           _one,
                           _constant,
                           _constant)

  rparams = RiskParameters{SEIR}(Float64[],
                                 [1.0],
                                 [2.5, 5.0],
                                 Float64[],
                                 [0.1],
                                 [0.05])

  sim = Simulation(pop,
                   [State_I; fill(State_S, n-1)],
                   rf,
                   rparams,
                   JC69(0.05))

  @test_nowarn simulate!(rng, sim, tmax=100.0)

  state_counts = [Pathogen._count_by_state(sim.events, State_S, 50.0)
                  Pathogen._count_by_state(sim.events, State_E, 50.0)
                  Pathogen._count_by_state(sim.events, State_I, 50.0)
                  Pathogen._count_by_state(sim.events, State_R, 50.0)]

  @test sum(state_counts) == n
  @test length(sim.disease_states) == n
  @test size(sim.transmission_network.internal) == (n, n)

  obs = observe(rng, sim, Uniform(), Uniform(), seq_len=100) #seq length of 100
  @test length(obs.infection) == n

  rpriors = RiskPriors{SEIR}(UnivariateDistribution[],
                             [Uniform(0.0, 2.0)],
                             [Uniform(0.0, 4.0); Uniform(1.0, 8.0)],
                             UnivariateDistribution[],
                             [Uniform(0.0, 1.0)],
                             [Uniform(0.0, 1.0)])
  smpriors = UnivariateDistribution[Uniform(0.0, 0.05)]
  ee = EventExtents{SEIR}((2.0, 20.0), 2.0, 2.0)
  mcmc = MCMC(obs, ee, pop, [State_I; fill(State_S, n-1)], rf, rpriors, JC69abs, smpriors)
  start!(mcmc, 3, attempts=100)
  @test length(mcmc.markov_chains) == 3
  @test all([length(mcmc.markov_chains[i].risk_parameters[1]) for i=1:3] .== length(rparams))

  # Check bounds of Events initialization
  for j = 1:3, i = 1:n
    if mcmc.markov_chains[j].events[1].exposure[i] > -Inf
      l, u = Pathogen._bounds(i, State_E, ee, obs, mcmc.markov_chains[j].events[1], mcmc.markov_chains[j].transmission_network[1])
      @test l < u
    end
    if mcmc.markov_chains[j].events[1].infection[i] > -Inf
      l, u = Pathogen._bounds(i, State_I, ee, obs, mcmc.markov_chains[j].events[1], mcmc.markov_chains[j].transmission_network[1])
      @test l < u
    end
    if mcmc.markov_chains[j].events[1].removal[i] > -Inf
      l, u = Pathogen._bounds(i, State_R, ee, obs, mcmc.markov_chains[j].events[1], mcmc.markov_chains[j].transmission_network[1])
      @test l < u
    end
  end
  iterate!(mcmc, 100, 1.0)
  @test all([mcmc.markov_chains[i].iterations for i = 1:3] .== 100)
  tnd = TNDistribution(mcmc, burnin=2)
  @test sum(tnd) ≈ sum(obs.infection .!== NaN)
  @test sum(mode(tnd)) ≈ sum(obs.infection .!== NaN)
  @test all(Pathogen._outdegree(tnd) .>= 0.0)
  @test typeof(summary(mcmc, burnin=50, thin=2)) == DataFrame

  tn = mcmc.markov_chains[1].transmission_network[end]
  for i = 1:n
    pt = Pathogen._pathway_to(i, tn, depth=1)
    if pt != nothing
      @test i in Pathogen._pathway_from(pt[2], tn, depth=1)
    end
  end
  @test minimum(mcmc.markov_chains[1].events[end]) > -Inf
  phylotree = generate(Tree, mcmc.markov_chains[1].events[end], obs, mcmc.markov_chains[1].transmission_network[end])[1]
  @test typeof(phylotree) == Tree
  _treeplot(phylotree)
  @test typeof(summary(mcmc)) == DataFrame
end


@testset "SEI PhyloILM" begin
  # Some commonly used functions/examples provided in helpers/RiskFunctions.jl
  # For SEIR, risk functions and parameters in order of: sparks, susceptibility, infectivity, transmissibility, latency, and removal
  rf = RiskFunctions{SEI}(_zero,
                           _coefficient,
                           _powerlaw,
                           _one,
                           _constant)

  rparams = RiskParameters{SEI}(Float64[],
                                 [1.0],
                                 [2.5, 5.0],
                                 Float64[],
                                 [0.1])

  sim = Simulation(pop,
                   [State_I; fill(State_S, n-1)],
                   rf,
                   rparams,
                   JC69(0.05))

  @test_nowarn simulate!(rng, sim, tmax=100.0)

  state_counts = [Pathogen._count_by_state(sim.events, State_S, 50.0)
                  Pathogen._count_by_state(sim.events, State_E, 50.0)
                  Pathogen._count_by_state(sim.events, State_I, 50.0)]

  @test sum(state_counts) == n
  @test length(sim.disease_states) == n
  @test size(sim.transmission_network.internal) == (n, n)

  obs = observe(rng, sim, Uniform(), seq_len=100) #seq length of 100
  @test length(obs.infection) == n

  rpriors = RiskPriors{SEI}(UnivariateDistribution[],
                             [Uniform(0.0, 2.0)],
                             [Uniform(0.0, 4.0); Uniform(1.0, 8.0)],
                             UnivariateDistribution[],
                             [Uniform(0.0, 1.0)])
  smpriors = UnivariateDistribution[Uniform(0.0, 0.05)]
  ee = EventExtents{SEI}((2.0, 20.0), 2.0)
  mcmc = MCMC(obs, ee, pop, [State_I; fill(State_S, n-1)], rf, rpriors, JC69abs, smpriors)
  start!(mcmc, 3, attempts=100)
  @test length(mcmc.markov_chains) == 3
  @test all([length(mcmc.markov_chains[i].risk_parameters[1]) for i=1:3] .== length(rparams))

  # Check bounds of Events initialization
  for j = 1:3, i = 1:n
    if mcmc.markov_chains[j].events[1].exposure[i] > -Inf
      l, u = Pathogen._bounds(i, State_E, ee, obs, mcmc.markov_chains[j].events[1], mcmc.markov_chains[j].transmission_network[1])
      @test l < u
    end
    if mcmc.markov_chains[j].events[1].infection[i] > -Inf
      l, u = Pathogen._bounds(i, State_I, ee, obs, mcmc.markov_chains[j].events[1], mcmc.markov_chains[j].transmission_network[1])
      @test l < u
    end
  end
  iterate!(mcmc, 100, 1.0)
  @test all([mcmc.markov_chains[i].iterations for i = 1:3] .== 100)
  tnd = TNDistribution(mcmc, burnin=2)
  @test sum(tnd) ≈ sum(obs.infection .!== NaN)
  @test sum(mode(tnd)) ≈ sum(obs.infection .!== NaN)
  @test all(Pathogen._outdegree(tnd) .>= 0.0)
  @test typeof(summary(mcmc, burnin=50, thin=2)) == DataFrame

  tn = mcmc.markov_chains[1].transmission_network[end]
  for i = 1:n
    pt = Pathogen._pathway_to(i, tn, depth=1)
    if pt != nothing
      @test i in Pathogen._pathway_from(pt[2], tn, depth=1)
    end
  end
  @test minimum(mcmc.markov_chains[1].events[end]) > -Inf
  phylotree = generate(Tree, mcmc.markov_chains[1].events[end], obs, mcmc.markov_chains[1].transmission_network[end])[1]
  @test typeof(phylotree) == Tree
  _treeplot(phylotree)
  @test typeof(summary(mcmc)) == DataFrame
end


@testset "SIR PhyloILM" begin
  # Some commonly used functions/examples provided in helpers/RiskFunctions.jl
  # For SEIR, risk functions and parameters in order of: sparks, susceptibility, infectivity, transmissibility, latency, and removal
  rf = RiskFunctions{SIR}(_zero,
                           _coefficient,
                           _powerlaw,
                           _one,
                           _constant)

  rparams = RiskParameters{SIR}(Float64[],
                                 [1.0],
                                 [2.5, 5.0],
                                 Float64[],
                                 [0.05])

  sim = Simulation(pop,
                   [State_I; fill(State_S, n-1)],
                   rf,
                   rparams,
                   JC69(0.05))

  @test_nowarn simulate!(rng, sim, tmax=100.0)

  state_counts = [Pathogen._count_by_state(sim.events, State_S, 50.0)
                  Pathogen._count_by_state(sim.events, State_I, 50.0)
                  Pathogen._count_by_state(sim.events, State_R, 50.0)]

  @test sum(state_counts) == n
  @test length(sim.disease_states) == n
  @test size(sim.transmission_network.internal) == (n, n)

  obs = observe(rng, sim, Uniform(), Uniform(), seq_len=100) #seq length of 100
  @test length(obs.infection) == n

  rpriors = RiskPriors{SIR}(UnivariateDistribution[],
                             [Uniform(0.0, 2.0)],
                             [Uniform(0.0, 4.0); Uniform(1.0, 8.0)],
                             UnivariateDistribution[],
                             [Uniform(0.0, 1.0)])
  smpriors = UnivariateDistribution[Uniform(0.0, 0.05)]
  ee = EventExtents{SIR}(2.0, 2.0)
  mcmc = MCMC(obs, ee, pop, [State_I; fill(State_S, n-1)], rf, rpriors, JC69abs, smpriors)
  start!(mcmc, 3, attempts=100)
  @test length(mcmc.markov_chains) == 3
  @test all([length(mcmc.markov_chains[i].risk_parameters[1]) for i=1:3] .== length(rparams))

  # Check bounds of Events initialization
  for j = 1:3, i = 1:n
    if mcmc.markov_chains[j].events[1].infection[i] > -Inf
      l, u = Pathogen._bounds(i, State_I, ee, obs, mcmc.markov_chains[j].events[1], mcmc.markov_chains[j].transmission_network[1])
      @test l < u
    end
    if mcmc.markov_chains[j].events[1].removal[i] > -Inf
      l, u = Pathogen._bounds(i, State_R, ee, obs, mcmc.markov_chains[j].events[1], mcmc.markov_chains[j].transmission_network[1])
      @test l < u
    end
  end
  iterate!(mcmc, 100, 1.0)
  @test all([mcmc.markov_chains[i].iterations for i = 1:3] .== 100)
  tnd = TNDistribution(mcmc, burnin=2)
  @test sum(tnd) ≈ sum(obs.infection .!== NaN)
  @test sum(mode(tnd)) ≈ sum(obs.infection .!== NaN)
  @test all(Pathogen._outdegree(tnd) .>= 0.0)
  @test typeof(summary(mcmc, burnin=50, thin=2)) == DataFrame

  tn = mcmc.markov_chains[1].transmission_network[end]
  for i = 1:n
    pt = Pathogen._pathway_to(i, tn, depth=1)
    if pt != nothing
      @test i in Pathogen._pathway_from(pt[2], tn, depth=1)
    end
  end
  @test minimum(mcmc.markov_chains[1].events[end]) > -Inf
  phylotree = generate(Tree, mcmc.markov_chains[1].events[end], obs, mcmc.markov_chains[1].transmission_network[end])[1]
  @test typeof(phylotree) == Tree
  _treeplot(phylotree)
  @test typeof(summary(mcmc)) == DataFrame
end


@testset "SI PhyloILM" begin
  # Some commonly used functions/examples provided in helpers/RiskFunctions.jl
  # For SEIR, risk functions and parameters in order of: sparks, susceptibility, infectivity, transmissibility, latency, and removal
  rf = RiskFunctions{SI}(_zero,
                           _coefficient,
                           _powerlaw,
                           _one)

  rparams = RiskParameters{SI}(Float64[],
                                 [1.0],
                                 [2.5, 5.0],
                                 Float64[])

  sim = Simulation(pop,
                   [State_I; fill(State_S, n-1)],
                   rf,
                   rparams,
                   JC69(0.05))

  @test_nowarn simulate!(rng, sim, tmax=100.0)

  state_counts = [Pathogen._count_by_state(sim.events, State_S, 50.0)
                  Pathogen._count_by_state(sim.events, State_I, 50.0)]

  @test sum(state_counts) == n
  @test length(sim.disease_states) == n
  @test size(sim.transmission_network.internal) == (n, n)

  obs = observe(rng, sim, Uniform(), seq_len=100) #seq length of 100
  @test length(obs.infection) == n

  rpriors = RiskPriors{SI}(UnivariateDistribution[],
                             [Uniform(0.0, 2.0)],
                             [Uniform(0.0, 4.0); Uniform(1.0, 8.0)],
                             UnivariateDistribution[])
  smpriors = UnivariateDistribution[Uniform(0.0, 0.05)]
  ee = EventExtents{SI}(2.0)
  mcmc = MCMC(obs, ee, pop, [State_I; fill(State_S, n-1)], rf, rpriors, JC69abs, smpriors)
  start!(mcmc, 3, attempts=100)
  @test length(mcmc.markov_chains) == 3
  @test all([length(mcmc.markov_chains[i].risk_parameters[1]) for i=1:3] .== length(rparams))

  # Check bounds of Events initialization
  for j = 1:3, i = 1:n
    if mcmc.markov_chains[j].events[1].infection[i] > -Inf
      l, u = Pathogen._bounds(i, State_I, ee, obs, mcmc.markov_chains[j].events[1], mcmc.markov_chains[j].transmission_network[1])
      @test l < u
    end
  end
  iterate!(mcmc, 100, 1.0)
  @test all([mcmc.markov_chains[i].iterations for i = 1:3] .== 100)
  tnd = TNDistribution(mcmc, burnin=2)
  @test sum(tnd) ≈ sum(obs.infection .!== NaN)
  @test sum(mode(tnd)) ≈ sum(obs.infection .!== NaN)
  @test all(Pathogen._outdegree(tnd) .>= 0.0)
  @test typeof(summary(mcmc, burnin=50, thin=2)) == DataFrame

  tn = mcmc.markov_chains[1].transmission_network[end]
  for i = 1:n
    pt = Pathogen._pathway_to(i, tn, depth=1)
    if pt != nothing
      @test i in Pathogen._pathway_from(pt[2], tn, depth=1)
    end
  end
  @test minimum(mcmc.markov_chains[1].events[end]) > -Inf
  phylotree = generate(Tree, mcmc.markov_chains[1].events[end], obs, mcmc.markov_chains[1].transmission_network[end])[1]
  @test typeof(phylotree) == Tree
  _treeplot(phylotree)
  @test typeof(summary(mcmc)) == DataFrame
end


# Shared fixtures use the risk functions loaded on every worker by runtests.jl.
function rng_fixture(::Type{S}, phylo) where S <: DiseaseStateSequence
  n = 5
  pop = Population(DataFrame(riskfactor1=ones(n)), ones(n, n))
  funcs = Any[_zero, _constant, _powerlaw, _one]
  params = [Float64[], [1.0], [2.0, 1.0], Float64[]]
  priors = [UnivariateDistribution[], UnivariateDistribution[Uniform(0.5, 1.5)],
            UnivariateDistribution[Uniform(1.0, 3.0), Uniform(0.5, 1.5)], UnivariateDistribution[]]
  for state in (State_E, State_R)
    if state ∈ S
      push!(funcs, _constant)
      push!(params, [0.5])
      push!(priors, UnivariateDistribution[Uniform(0.3, 0.7)])
    end
  end
  rf = RiskFunctions{S}(funcs...)
  rp = RiskParameters{S}(params...)
  states = [State_I; fill(State_S, n-1)]
  sim = phylo ? Simulation(pop, states, rf, rp, JC69(0.05)) : Simulation(pop, states, rf, rp)
  extents = S == SEIR ? EventExtents{S}(30.0, 0.2, 0.2) :
            S == SEI ? EventExtents{S}(30.0, 0.2) :
            S == SIR ? EventExtents{S}(0.2, 0.2) : EventExtents{S}(0.2)
  return sim, RiskPriors{S}(priors...), extents, states
end

function rng_observe(rng, sim::Simulation{S}) where S
  if State_R ∈ S
    observe(rng, sim, Uniform(0.0, 0.1), Uniform(0.0, 0.1); seq_len=16, force=true)
  else
    observe(rng, sim, Uniform(0.0, 0.1); seq_len=16)
  end
end

function rng_mcmc(obs, sim, priors, extents, states)
  if sim.substitution_model === nothing
    MCMC(obs, extents, sim.population, states, sim.risk_functions, priors)
  else
    MCMC(obs, extents, sim.population, states, sim.risk_functions, priors,
         JC69abs, UnivariateDistribution[Uniform(0.01, 0.1)])
  end
end

function rng_snapshot(mc)
  (mc.iterations, [convert(Array{Float64, 2}, e) for e in mc.events],
   [convert(Vector{Float64}, rp) for rp in mc.risk_parameters],
   [(tn.external, tn.internal) for tn in mc.transmission_network],
   mc.substitution_model === nothing ? nothing :
     [Tuple(getproperty(sm, f) for f in propertynames(sm)) for sm in mc.substitution_model],
   mc.log_posterior)
end

@testset "Explicit RNGs" begin
  generators = (MersenneTwister, StableRNG)
  for S in (SEIR, SEI, SIR, SI), phylo in (false, true), R in generators
    @testset "$S phylo=$phylo $R" begin
      sim, priors, extents, states = rng_fixture(S, phylo)
      # Every simulation and observation draw must come from the supplied RNG.
      Random.seed!(908)
      expected_default = rand()
      Random.seed!(908)
      a, b = deepcopy(sim), deepcopy(sim)
      ra, rb = R(42), R(42)
      simulate!(ra, a; tmax=50.0)
      simulate!(rb, b; tmax=50.0)
      @test isequal(convert(Array{Float64, 2}, a.events), convert(Array{Float64, 2}, b.events))
      oa, ob = rng_observe(ra, a), rng_observe(rb, b)
      @test isequal(oa.infection, ob.infection)
      @test isequal(oa.removal, ob.removal)
      @test isequal(oa.seq, ob.seq)
      @test rand() == expected_default
      @test rand(ra) == rand(rb)
      @test rand(copy(ra)) != rand(R(42))

      ma = rng_mcmc(oa, a, priors, extents, states)
      mb = deepcopy(ma)
      # Exercise one-chain delegation, initialization, and MCMC proposals.
      Random.seed!(908)
      start!(ra, ma, 1; attempts=20)
      start!(rb, mb; attempts=20)
      @test length(ma.markov_chains) == length(mb.markov_chains) == 1
      @test isequal(rng_snapshot(ma.markov_chains[1]), rng_snapshot(mb.markov_chains[1]))
      iterate!(ra, ma, 4, 0.1; adapt_cov=0)
      iterate!(rb, mb, 4, 0.1; adapt_cov=0)
      @test isequal(rng_snapshot(ma.markov_chains[1]), rng_snapshot(mb.markov_chains[1]))
      @test rand() == expected_default
      @test rand(ra) == rand(rb)
      # Cover conditional event proposals and coupled network augmentation.
      opts = phylo ? (DA_couple=true,) : (condition_on_network=true,)
      iterate!(ra, ma.markov_chains[1], ma, 2, 0.1; adapt_cov=0, opts...)
      iterate!(rb, mb.markov_chains[1], mb, 2, 0.1; adapt_cov=0, opts...)
      @test isequal(rng_snapshot(ma.markov_chains[1]), rng_snapshot(mb.markov_chains[1]))

      # The default-RNG wrappers retain the same serial semantics.
      mc = deepcopy(ma.markov_chains[1])
      md = deepcopy(mc)
      Random.seed!(123)
      iterate!(mc, ma, 1, 0.1; adapt_cov=0)
      Random.seed!(123)
      iterate!(Pathogen._default_rng(), md, ma, 1, 0.1; adapt_cov=0)
      @test isequal(rng_snapshot(mc), rng_snapshot(md))
    end
  end

  @testset "Distributed streams" for phylo in (false, true)
    sim, priors, extents, states = rng_fixture(SIR, phylo)
    simulate!(MersenneTwister(42), sim; tmax=50.0)
    obs = rng_observe(MersenneTwister(13), sim)
    base = rng_mcmc(obs, sim, priors, extents, states)
    a, b = deepcopy(base), deepcopy(base)
    ra, rb = MersenneTwister(77), MersenneTwister(77)
    Random.seed!(908)
    expected_default = rand()
    Random.seed!(908)
    start!(ra, a, 2; attempts=20)
    start!(rb, b, 2; attempts=20)
    @test length(a.markov_chains) == length(b.markov_chains) == 2
    @test isequal(rng_snapshot.(a.markov_chains), rng_snapshot.(b.markov_chains))
    @test !isequal(rng_snapshot(a.markov_chains[1]), rng_snapshot(a.markov_chains[2]))
    iterate!(ra, a, 4, 0.1; adapt_cov=0)
    iterate!(rb, b, 4, 0.1; adapt_cov=0)
    @test isequal(rng_snapshot.(a.markov_chains), rng_snapshot.(b.markov_chains))
    @test rand(ra) == rand(rb)
    @test rand() == expected_default
  end
end


include("inference_performance.jl")

include("risk_specialization.jl")

include("event_proposals.jl")
