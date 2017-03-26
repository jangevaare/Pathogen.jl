__precompile__()

module Pathogen

  # Dependencies
  using DataFrames
  using PhyloTrees
  using PhyloModels
  using Distributions
  using RecipesBase
  using ProgressMeter

  # Functions to be extended
  import
    Base.getindex,
    Base.push!,
    Base.append!,
    Base.rand,
    Base.length,
    Base.convert,
    Base.size,
    Base.show,
    Base.Array,
    Base.copy,
    Base.deleteat!,
    PhyloModels.simulate!,
    PhyloModels.loglikelihood,
    PhyloModels.logprior,
    PhyloModels.propose

  # Source files
  ## Core
  include("core/risks.jl")
  include("core/states.jl")
  include("core/rates.jl")
  include("core/events.jl")
  include("core/networks.jl")
  include("core/observe.jl")

  ## Utilities
  include("utilities/pathways.jl")
  include("utilities/states.jl")
  include("utilities/trees.jl")
  include("utilities/plotting.jl")

  ## Simulation
  include("simulation.jl")

  ## Inference
  include("inference/risks.jl")
  include("inference/events.jl")
  include("inference/networks.jl")
  include("inference/loglikelihood.jl")
  include("inference/mcmc.jl")

  ## Visualization
  include("plotrecipes.jl")

  # New types and functions
  export
    ## Core
    RiskFunctions,
      SEIR_RiskFunctions,
      SIR_RiskFunctions,
      SEI_RiskFunctions,
      SI_RiskFunctions,
    RiskParameters,
      SEIR_RiskParameters,
      SIR_RiskParameters,
      SEI_RiskParameters,
      SI_RiskParameters,
    NetworkRates,
    Rates,
      SEIR_Rates,
      SIR_Rates,
      SEI_Rates,
      SI_Rates,
    Network,
    States,
      SEIR_States,
      SIR_States,
      SEI_States,
      SI_States,
    Events,
      SEIR_Events,
      SIR_Events,
      SEI_Events,
      SI_Events,
    EventObservations,
      SEIR_EventObservations,
      SIR_EventObservations,
      SEI_EventObservations,
      SI_EventObservations,

    ## Utilities
    generate_tree,
    pathwayto,
    pathwayfrom,

    ## Simulation
    initialize_simulation,
    simulate!,
    observe,

    ## Inference
    RiskParameterPriors,
    generate_events,
    PathogenTrace,
    PathogenIteration,
    initialize_mcmc,
    mcmc!,
    propose
end
