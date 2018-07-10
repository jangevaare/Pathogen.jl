using Base.Test
using Pathogen
using DataFrames
using Distributions

include(Pkg.dir("Pathogen")*"/examples/risk_functions.jl")

# Set RNG seed
srand(5432)

# Define population
n = 25
x_coordinates = rand(Uniform(0, 5), n)
y_coordinates = rand(Uniform(0, 5), n)
riskfactor1 = rand(Gamma(), n)
pop = DataFrame(x = x_coordinates,
                y = y_coordinates,
                riskfactor1 = riskfactor1)


@testset "SEIR Model" begin
  # Some commonly used functions/examples provided in helpers/RiskFunctions.jl
  # For SEIR, risk functions and parameters in order of: sparks, susceptibility, transmissibility, infectivity, latency, and removal
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

  sim = Simulation(pop, rf, rparams)

  simulate!(sim, tmax=100.0)

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
  mcmc = MCMC(obs, ee, pop, rf, rpriors)
  start!(mcmc, markov_chains=3)
  @test length(mcmc.markov_chains) == 3
  @test all([length(mcmc.markov_chains[i].risk_parameters[1]) for i=1:3] .== length(rparams))
  while all([mcmc.markov_chains[i].iterations for i=1:3] .< 100)
    next!(mcmc, diagm([0.0001; 0.01; 0.01; 0.01; 0.001; 0.001]), 0.5)
  end
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

  sim = Simulation(pop, rf, rparams)

  simulate!(sim, tmax=100.0)

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
  start!(mcmc, markov_chains=3)
  @test length(mcmc.markov_chains) == 3
  @test all([length(mcmc.markov_chains[i].risk_parameters[1]) for i=1:3] .== length(rparams))
  while all([mcmc.markov_chains[i].iterations for i=1:3] .< 100)
    next!(mcmc, diagm([0.0001; 0.01; 0.01; 0.01; 0.001]), 0.5)
  end
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

  sim = Simulation(pop, rf, rparams)

  simulate!(sim, tmax=100.0)

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
  start!(mcmc, markov_chains=3)
  @test length(mcmc.markov_chains) == 3
  @test all([length(mcmc.markov_chains[i].risk_parameters[1]) for i=1:3] .== length(rparams))
  while all([mcmc.markov_chains[i].iterations for i=1:3] .< 100)
    next!(mcmc, diagm([0.0001; 0.01; 0.01; 0.01; 0.001]), 0.5)
  end
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

  sim = Simulation(pop, rf, rparams)

  simulate!(sim, tmax=100.0)

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
  start!(mcmc, markov_chains=3)
  @test length(mcmc.markov_chains) == 3
  @test all([length(mcmc.markov_chains[i].risk_parameters[1]) for i=1:3] .== length(rparams))
  while all([mcmc.markov_chains[i].iterations for i=1:3] .< 100)
    next!(mcmc, diagm([0.0001; 0.01; 0.01; 0.01]), 0.5)
  end
end
