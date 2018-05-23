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


@testset "SEIR Simulation" begin
  # Some commonly used functions/examples provided in helpers/RiskFunctions.jl
  # For SEIR, risk functions and parameters in order of: sparks, susceptibility, transmissibility, infectivity, latency, and removal
  rf = RiskFunctions{SEIR}(Pathogen._constant,
                           Pathogen._coefficient,
                           Pathogen._powerlaw,
                           Pathogen._one,
                           Pathogen._constant,
                           Pathogen._constant)

  rp = RiskParameters{SEIR}([0.001],
                            [1.0],
                            [2.0, 5.0],
                            Float64[],
                            [0.1],
                            [0.05])

  sim = Simulation(pop, rf, rp)

  simulate!(sim, nmax=50)

  @test length(sim.disease_states) == n
  @test size(sim.transmission_network.internal) == (n, n)
end

@testset "SEI Simulation" begin
  # Some commonly used functions/examples provided in helpers/RiskFunctions.jl
  # For SEIR, risk functions and parameters in order of: sparks, susceptibility, transmissibility, infectivity, latency, and removal
  rf = RiskFunctions{SEI}(Pathogen._constant,
                          Pathogen._coefficient,
                          Pathogen._powerlaw,
                          Pathogen._one,
                          Pathogen._constant)

  rp = RiskParameters{SEI}([0.001],
                           [1.0],
                           [2.0, 5.0],
                           Float64[],
                           [0.1])

  sim = Simulation(pop, rf, rp)

  simulate!(sim, nmax=50)

  @test length(sim.disease_states) == n
  @test size(sim.transmission_network.internal) == (n, n)
end

@testset "SIR Simulation" begin
  # Some commonly used functions/examples provided in helpers/RiskFunctions.jl
  # For SEIR, risk functions and parameters in order of: sparks, susceptibility, transmissibility, infectivity, latency, and removal
  rf = RiskFunctions{SIR}(Pathogen._constant,
                           Pathogen._coefficient,
                           Pathogen._powerlaw,
                           Pathogen._one,
                           Pathogen._constant)

  rp = RiskParameters{SIR}([0.001],
                           [1.0],
                           [2.0, 5.0],
                           Float64[],
                           [0.05])

  sim = Simulation(pop, rf, rp)

  simulate!(sim, nmax=50)

  @test length(sim.disease_states) == n
  @test size(sim.transmission_network.internal) == (n, n)
end

@testset "SI Simulation" begin
  # Some commonly used functions/examples provided in helpers/RiskFunctions.jl
  # For SEIR, risk functions and parameters in order of: sparks, susceptibility, transmissibility, infectivity, latency, and removal
  rf = RiskFunctions{SI}(Pathogen._constant,
                         Pathogen._coefficient,
                         Pathogen._powerlaw,
                         Pathogen._one)

  rp = RiskParameters{SI}([0.001],
                          [1.0],
                          [2.0, 5.0],
                          Float64[])

  sim = Simulation(pop, rf, rp)

  simulate!(sim, nmax=50, pmax=60.0)

  @test length(sim.disease_states) == n
  @test size(sim.transmission_network.internal) == (n, n)
end
