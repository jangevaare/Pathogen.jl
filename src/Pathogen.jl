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
  findstate,
  plotdata,
  convert,

  # substitute.jl
  jc69,

  # simulate.jl
  create_seq,
  create_population,
  create_powerlaw,
  create_constantrate,
  create_ratearray,
  onestep!,

  # infer.jl
  SEIR_surveilance
  create_tree,
  seqdistance,
  branchloglikelihood,
  create_logprior1

include("types.jl")
include("utilities.jl")
include("substitute.jl")
include("simulate.jl")
include("infer.jl")

end
