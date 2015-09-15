module Pathogen

using Distributions, Distances, DataFrames, BioSeq, ProgressMeter

export
  # types.jl
  Population,
  RateArray,
  SEIR_actual,
  SEIR_observed,
  SEIR_augmented,

  Priors,
    ILM_priors,
      SEIR_priors,
    Detection_priors,
      Lag_priors,
    Mutation_priors,
      JC69_priors,

  Trace,
    ILM_trace,
      SEIR_trace,
    Detection_trace,
      Lag_trace,
    Mutation_trace,
      JC69_trace,

  # utilities.jl
  findstate,
  plotdata,
  convert,
  maximum,
  isseq,

  # mutate.jl
  jc69,

  # simulate.jl
  create_seq,
  create_population,
  create_powerlaw,
  create_constantrate,
  create_ratearray,
  onestep!,

  # infer.jl
  surveil,
  augment,
  logprior,
  rand_prior,
  seq_distances,
  network_loglikelihood,
  SEIR_loglikelihood,
  initialize,
  MCMC

include("types.jl")
include("utilities.jl")
include("mutate.jl")
include("simulate.jl")
include("infer.jl")

end
