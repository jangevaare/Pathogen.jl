module Pathogen


using Distributions, Distances, DataFrames, BioSeq, ProgressMeter, UnicodePlots


export
  # types.jl
  Population,

  RateArray,

  Actual,
    SEIR_actual,
    SIR_actual,

  Observed,
    SEIR_observed,
    SIR_observed,

  Augmented,
    SEIR_augmented,
    SIR_augmented,

  Priors,
    ILM_priors,
      SEIR_priors,
      SIR_priors,
    Detection_priors,
      Lag_priors,
    Mutation_priors,
      JC69_priors,

  Trace,
    ILM_trace,
      SEIR_trace,
      SIR_trace,
    Detection_trace,
      Lag_trace,
    Mutation_trace,
      JC69_trace,

  # utilities.jl
  findstate,
  findnetwork,
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

  # inference
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
  SIR_loglikelihood,
  initialize,
  MCMC


include("types.jl")
include("utilities.jl")
include("mutate.jl")
include("simulate.jl")
include("inference/eventtimes/seir.jl")
include("inference/eventtimes/sir.jl")
include("inference/ilm.jl")
include("inference/network.jl")
include("inference/mcmc/common.jl")
include("inference/mcmc/seir.jl")
include("inference/mcmc/sir.jl")


end
