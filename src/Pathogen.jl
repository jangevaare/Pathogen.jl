module Pathogen

  # Dependencies
  using DataFrames
  using PhyloTrees
  using Distributions
  using RecipesBase
  using ProgressMeter

  # Functions to be extended
  import
    Base.push!,
    Base.append!,
    Base.rand,
    Base.length,
    Base.Array,
    Base.Vector,
    Base.size,
    Base.copy

  # Source files
  ## Core
  include("core/risks.jl")
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
  include("simulation/initialization.jl")
  include("simulation/run.jl")

  ## Inference
  include("inference/riskparameters.jl")
  include("inference/eventtimes.jl")
  include("inference/exposurenetwork.jl")
  include("inference/loglikelihood.jl")
  include("inference/mcmc.jl")

  ## Visualization
  include("plotrecipes.jl")

  # New types and functions
  export
    ## Core
    RiskFunctions,
    RiskParameters,
    Rates,
    Events,
    Network,

    ## Utilities
    pathwayto,
    pathwayfrom,
    findstate,
    generatetree,

    ## Simulation
    initialize_rates,
    initialize_simulation,
    update_events!,
    update_rates!,
    simulate!,
    EventObservations,
    observe,

    ## Inference
    RiskParameterPriors,
    EventPriors,
    logprior,
    PathogenTrace,
    PathogenIteration,
    PathogenProposal,
    loglikelihood
end
