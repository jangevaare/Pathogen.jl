module Pathogen

using Distances, Distributions, DataFrames, BioSeq, ProgressMeter

export
  # types.jl
  population,
  Tree,
  TreeNode,
  TreeLeaf,
  TreeFeature,

  # utilities.jl
  generate_seq,
  generate_2bitseq,
  findstate,
  plotdata,
  convert,
  surveil,

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
