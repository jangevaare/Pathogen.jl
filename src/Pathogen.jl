module Pathogen


using Distributions, Distances, DataFrames, BioSeq, ProgressMeter, UnicodePlots


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
    Posterior_trace

  # utilities.jl
  findstate,
  plotdata,
  convert,
  maximum,
  isseq,
  pathwayto,
  pathwaysto,
  pathwayfrom,
  pathwaysfrom,

  # mutate.jl
  jc69q,
  jc69p,

  # simulate.jl
  create_seq,
  create_population,
  create_powerlaw,
  create_constantrate,
  create_ratearray,
  onestep!,

  # infer.jl
  surveil,
  propose_augment,
  logprior,
  rand_prior,
  propose_network,
  seq_distances,
  phylogenetic_network_loglikelihood,
  exposure_network_loglikelihood,
  detection_loglikelihood,
  SEIR_loglikelihood,
  initialize,
  MCMC


include("types.jl")
include("utilities.jl")
include("mutate.jl")
include("simulate.jl")
include("infer.jl")


end
