module Pathogen

using Distances, Distributions, DataFrames, BioSeq, ProgressMeter

export
  # types.jl
  population,

  # utilities.jl
  GenerateSequence,

  # substitute.jl
  JC69,

  # simulate.jl
  CreatePopulation,
  CreatePowerLaw,
  CreateConstantRate,
  CreateRateArray

include("types.jl")
include("utilities.jl")
include("substitute.jl")
include("simulate.jl")

end
