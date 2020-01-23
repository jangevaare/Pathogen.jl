using Test, Distributed, Random, LinearAlgebra, Distances, DataFrames, Distributions, Pathogen

addprocs(3)

@everywhere using Pathogen
@everywhere include(joinpath(@__DIR__, "risk_functions.jl"))

# Set RNG seed
Random.seed!(5432)

# using Logging
# global_logger(ConsoleLogger(stderr, LogLevel(-10000)))

# Define population
n = 25
risks = DataFrame(x = rand(Uniform(0, 5), n),
                  y = rand(Uniform(0, 5), n),
                  riskfactor1 = rand(Gamma(), n))

pop = Population(risks,
                 [euclidean([risks[i, :x]; risks[i, :y]], [risks[j, :x]; risks[j, :y]]) for i = 1:n, j = 1:n])

@testset "SEIR Model" begin
  # Some commonly used functions/examples provided in helpers/RiskFunctions.jl
  # For SEIR, risk functions and parameters in order of: sparks, susceptibility, infectivity, transmissibility, latency, and removal
  rf = RiskFunctions{SEIR}(_constant,
                           _coefficient,
                           _powerlaw,
                           _one,
                           _constant,
                           _constant)

  rparams = RiskParameters{SEIR}([0.001],
                                 [1.0],
                                 [2.0, 5.0],
                                 Float64[],
                                 [0.1],
                                 [0.05])

  # sim = Simulation(pop, rf, rparams)
  sim = Simulation(pop,
                   [State_I; fill(State_S, n-1)],
                   0.0,
                   rf,
                   rparams)

  simulate!(sim, tmax=100.0)

  state_counts = [Pathogen._count_by_state(sim.events, State_S, 50.0)
                  Pathogen._count_by_state(sim.events, State_E, 50.0)
                  Pathogen._count_by_state(sim.events, State_I, 50.0)
                  Pathogen._count_by_state(sim.events, State_R, 50.0)]

  @test sum(state_counts) == n
  @test length(sim.disease_states) == n
  @test size(sim.transmission_network.internal) == (n, n)

  obs = observe(sim, Uniform(), Uniform())
  @test length(obs.infection) == n

  rpriors = RiskPriors{SEIR}([Uniform(0.0, 0.01)],
                             [Uniform(0.0, 2.0)],
                             [Uniform(0.0, 4.0); Uniform(1.0, 8.0)],
                             UnivariateDistribution[],
                             [Uniform(0.0, 1.0)],
                             [Uniform(0.0, 1.0)])
  ee = EventExtents{SEIR}(20.0, 2.0, 2.0)
  mcmc = MCMC(obs, ee, pop, [State_I; fill(State_S, n-1)], rf, rpriors)
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
  iterate!(mcmc, 100, 0.5)
  @test all([mcmc.markov_chains[i].iterations for i = 1:3] .== 100)
  @test sum(TNPosterior(mcmc)) ≈ sum(obs.infection .< Inf)-1
end

@testset "SEI Model" begin
  rf = RiskFunctions{SEI}(_constant,
                          _coefficient,
                          _powerlaw,
                          _one,
                          _constant)

  rparams = RiskParameters{SEI}([0.001],
                                [1.0],
                                [2.0, 5.0],
                                Float64[],
                                [0.1])

  sim = Simulation(pop,
                   [State_I; fill(State_S, n-1)],
                   0.0,
                   rf,
                   rparams)

  simulate!(sim, tmax=100.0)

  state_counts = [Pathogen._count_by_state(sim.events, State_S, 50.0)
                  Pathogen._count_by_state(sim.events, State_E, 50.0)
                  Pathogen._count_by_state(sim.events, State_I, 50.0)]

  @test sum(state_counts) == n

  @test length(sim.disease_states) == n
  @test size(sim.transmission_network.internal) == (n, n)

  obs = observe(sim, Uniform())
  @test length(obs.infection) == n

  rpriors = RiskPriors{SEI}([Uniform(0.0, 0.01)],
                            [Uniform(0.0, 2.0)],
                            [Uniform(0.0, 4.0); Uniform(1.0, 8.0)],
                            UnivariateDistribution[],
                            [Uniform(0.0, 1.0)])

  ee = EventExtents{SEI}(20.0, 2.0)
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
  iterate!(mcmc, 100, 0.5)
  @test all([mcmc.markov_chains[i].iterations for i = 1:3] .== 100)
  @test sum(TNPosterior(mcmc)) ≈ sum(obs.infection .< Inf)-1
end

@testset "SIR Model" begin
  rf = RiskFunctions{SIR}(_constant,
                          _coefficient,
                          _powerlaw,
                          _one,
                          _constant)

  rparams = RiskParameters{SIR}([0.001],
                                [1.0],
                                [2.0, 5.0],
                                Float64[],
                                [0.05])

  sim = Simulation(pop,
                   [State_I; fill(State_S, n-1)],
                   0.0,
                   rf,
                   rparams)

  simulate!(sim, tmax=100.0)

  state_counts = [Pathogen._count_by_state(sim.events, State_S, 50.0)
                  Pathogen._count_by_state(sim.events, State_I, 50.0)
                  Pathogen._count_by_state(sim.events, State_R, 50.0)]

  @test sum(state_counts) == n

  @test length(sim.disease_states) == n
  @test size(sim.transmission_network.internal) == (n, n)

  obs = observe(sim, Uniform(), Uniform())
  @test length(obs.infection) == n

  rpriors = RiskPriors{SIR}([Uniform(0.0, 0.01)],
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
  iterate!(mcmc, 100, 0.5)
  @test all([mcmc.markov_chains[i].iterations for i = 1:3] .== 100)
  @test sum(TNPosterior(mcmc)) ≈ sum(obs.infection .< Inf)-1
end

@testset "SI Model" begin
  rf = RiskFunctions{SI}(_constant,
                         _coefficient,
                         _powerlaw,
                         _one)

  rparams = RiskParameters{SI}([0.001],
                               [1.0],
                               [2.0, 5.0],
                               Float64[])

  sim = Simulation(pop,
                   [State_I; fill(State_S, n-1)],
                   0.0,
                   rf,
                   rparams)

  simulate!(sim, tmax=100.0)

  state_counts = [Pathogen._count_by_state(sim.events, State_S, 50.0)
                  Pathogen._count_by_state(sim.events, State_I, 50.0)]

  @test sum(state_counts) == n

  @test length(sim.disease_states) == n
  @test size(sim.transmission_network.internal) == (n, n)

  obs = observe(sim, Uniform())
  @test length(obs.infection) == n

  rpriors = RiskPriors{SI}([Uniform(0.0, 0.01)],
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
  iterate!(mcmc, 100, 0.5)
  @test all([mcmc.markov_chains[i].iterations for i = 1:3] .== 100)
  @test sum(TNPosterior(mcmc)) ≈ sum(obs.infection .< Inf)-1
end
