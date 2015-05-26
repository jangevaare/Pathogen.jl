module Pathogen

using Distributions, DataFrames, BioSeq

export
  # types.jl
  population,

  # utilities.jl
  generate_sequence,

  # substitute.jl
  JC69,

  # simulate.jl
  create_population,
  rate_array

include("types.jl")
include("utilities.jl")
include("substitute.jl")
include("simulate.jl")

end
