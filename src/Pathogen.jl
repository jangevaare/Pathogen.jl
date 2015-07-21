module Pathogen

using PDMats, Distributions, Distances, DataFrames, BioSeq, ProgressMeter

export
  # types.jl
  Population,
  RateArray,
  Tree,
  TreeVertex,
  TreeEdge,
  SEIR_actual,
  SEIR_observed,
  SEIR_augmented,
  SEIR_priors,
  SEIR_trace,

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
  SEIR_surveilance,
  SEIR_augmentation,
  SEIR_loglikelihood,
  SEIR_initialize,
  SEIR_MCMC,
  create_tree,
  seqdistance,
  branchloglikelihood

include("types.jl")
include("utilities.jl")
include("substitute.jl")
include("simulate.jl")
include("infer.jl")

end
