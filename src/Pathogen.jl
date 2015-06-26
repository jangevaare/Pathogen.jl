module Pathogen

using Distances, Distributions, DataFrames, BioSeq, ProgressMeter

export
  # types.jl
  Population,
  RateArray,
  Tree,
  TreeFeature,
  TreeLeaf,
  TreeNode,
  Tree2,

  # utilities.jl
  generate_seq,
  geneticdistance,
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
  onestep!,

  # infer.jl
  branchloglikelihood,
  treedistance,
  seqdistance

include("types.jl")
include("utilities.jl")
include("substitute.jl")
include("simulate.jl")
include("infer.jl")

end
