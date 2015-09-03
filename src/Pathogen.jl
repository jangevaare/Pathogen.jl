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
    ILM_priors
      SEIR_priors,
    Detection_priors,
      Lag_priors,
    Mutation_priors,
      JC69_priors

  Trace,
    SEIR_trace,
      ILM_trace,
    Detection_trace,
      Lag_trace,
    Mutation_trace,
      JC69_trace,

  # utilities.jl
  findstate,
  plotdata,
  convert,
  maximum,

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
  ILM_logprior,
  mutation_logprior,
  detection_logprior,
  SEIR_loglikelihood,
  SEIR_logprior,
  SEIR_initialize,
  SEIR_MCMC,
  seq_distances,
  seq_loglikelihood,
  network_loglikelihood

include("types.jl")
include("utilities.jl")
include("mutate.jl")
include("simulate.jl")
include("infer.jl")

end
