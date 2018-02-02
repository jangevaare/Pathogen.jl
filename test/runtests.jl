using Base.Test
using Pathogen
using DataFrames
using Distributions
using PhyloTrees
using PhyloModels

include("population.jl")
include("riskfunctions.jl")
include("SEIR_simulation.jl")
include("SEIR_inference.jl")
include("SIR_simulation.jl")
include("SIR_inference.jl")
include("SEI_simulation.jl")
include("SEI_inference.jl")
include("SI_simulation.jl")
include("SI_inference.jl")
