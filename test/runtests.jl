using Base.Test
using Pathogen
using DataFrames
using Distributions

# Set RNG seed
srand(5432)

# Define population
x_coordinates = rand(Uniform(0, 5), 25)
y_coordinates = rand(Uniform(0, 5), 25)
riskfactor1 = rand(Gamma(), 25)
pop = DataFrame(x = x_coordinates,
                y = y_coordinates,
                riskfactor1 = riskfactor1)



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

test_sim = Simulation(pop, rf, rp)
