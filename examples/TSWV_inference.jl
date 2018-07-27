#addprocs(3)
@everywhere using Pathogen
@everywhere using DataFrames, Distributions
using Plots, CSV

srand(5432)

# Use CSV.jl for DataFrames I/O
#
# We know the types of the columns, so we'll manually specify those.
# * Individual IDs are `Int64`
# * X,Y coordinates are `Float64`s
pop = CSV.read(Pkg.dir("Pathogen")*"/examples/02_TSWV_locations.csv", types=[Int64; Float64; Float64])

# Use julia's included CSV interface for simple vector of observation times
raw_observations = readcsv(Pkg.dir("Pathogen")*"/examples/02_TSWV_infection_observations.csv")[:]

# Create an `EventObservations` object with `Pathogen.jl`
obs = EventObservations{SI}(raw_observations)

@everywhere include(Pkg.dir("Pathogen")*"/examples/risk_functions.jl")

rf = RiskFunctions{SI}(_zero, # sparks function - we will assume no exogenous transmissions and set this to zero
                       _one, # susceptibility function - we do not have individual level risk factor information to explore here, so will set to a constant 1
                       _powerlaw_w_intercept, # transmissability function - we will use a powerlaw (with intercept) kernel. This provides a spatial and non-spatial component to infection transmissions. This has 3 parameters.
                       _one) # infectivity function - we do not have individual level risk factor information to explore here, so will set to a constant 1

rpriors = RiskPriors{SI}(UnivariateDistribution[], # empty `UnivariateDistribution` vector for all parameter-less functions
                    UnivariateDistribution[],
                    [Gamma(10.0, 10.0); Gamma(10.0, 10.0); Gamma(1.0, 1.0)], # Relatively uninformative priors with appropriate support
                    UnivariateDistribution[])

ee = EventExtents{SI}(14.0)

mcmc = MCMC(obs, ee, pop, rf, rpriors)
start!(mcmc, markov_chains = 3, attempts = 5000)

iterate!(mcmc, 10000, diagm([0.1; 0.1; 0.01]), 1.0)
