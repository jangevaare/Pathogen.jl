using Base.Test
using Pathogen
using DataFrames
using Distributions

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
  rf = RiskFunctions{SEIR}(Pathogen._constant,
                           Pathogen._coefficient,
                           Pathogen._powerlaw,
                           Pathogen._one,
                           Pathogen._constant,
                           Pathogen._constant)

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
end

@testset "SEI Model" begin
  # Some commonly used functions/examples provided in helpers/RiskFunctions.jl
  # For SEIR, risk functions and parameters in order of: sparks, susceptibility, transmissibility, infectivity, latency, and removal

  rf = RiskFunctions{SEI}(Pathogen._constant,
                        Pathogen._coefficient,
                        Pathogen._powerlaw,
                        Pathogen._one,
                        Pathogen._constant)

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
end

@testset "SIR Model" begin
  # Some commonly used functions/examples provided in helpers/RiskFunctions.jl
  # For SEIR, risk functions and parameters in order of: sparks, susceptibility, transmissibility, infectivity, latency, and removal
  rf = RiskFunctions{SIR}(Pathogen._constant,
                           Pathogen._coefficient,
                           Pathogen._powerlaw,
                           Pathogen._one,
                           Pathogen._constant)

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
end

@testset "SI Model" begin
  # Some commonly used functions/examples provided in helpers/RiskFunctions.jl
  # For SEIR, risk functions and parameters in order of: sparks, susceptibility, transmissibility, infectivity, latency, and removal
  rf = RiskFunctions{SI}(Pathogen._constant,
                         Pathogen._coefficient,
                         Pathogen._powerlaw,
                         Pathogen._one)

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
end
