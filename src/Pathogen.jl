module Pathogen

using Distances, Distributions, DataFrames, BioSeq, ProgressMeter

export
  # types.jl
  population,

  # utilities.jl
  generate_sequence,
  findstate,
  plotdata,

  # substitute.jl
  jc69,

  # simulate.jl
  create_population,
  create_powerlaw,
  create_constantrate,
  create_ratearray,
  onestep!

include("types.jl")
include("utilities.jl")
include("substitute.jl")
include("simulate.jl")

end
