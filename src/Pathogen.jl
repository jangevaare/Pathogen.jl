module Pathogen

using Distances, Graphs, Distributions, DataFrames, BioSeq, ProgressMeter

export
  # types.jl
  Population,
  RateArray,
  Tree,
  TreeVertex,
  TreeEdge,

  # utilities.jl
  generate_seq,
  seqdistance,
  findstate,
  plotdata,
  convert,
  seqsurveil,
  simplesurveil,
  create_tree,

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
  create_logprior1

include("types.jl")
include("utilities.jl")
include("substitute.jl")
include("simulate.jl")
include("infer.jl")

end
