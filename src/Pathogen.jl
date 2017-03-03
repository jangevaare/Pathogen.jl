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
    RiskParameters,
    Network,
    NetworkRates,
    Events,

    ## Utilities
    generate_tree,
    pathwayto,
    pathwayfrom,

    ## Simulation
    initialize_simulation,
    simulate!,
    EventObservations,
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
