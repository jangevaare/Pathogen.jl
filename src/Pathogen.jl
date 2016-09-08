module Pathogen

  # Dependencies
  using DataFrames
  using PhyloTrees
  using Distributions
  using RecipesBase

  # Functions to be extended
  import
    Base.push!,
    Base.append!,
    Base.rand

  # Source files
  ## Core
  include("core/risks.jl")
  include("core/rates.jl")
  include("core/events.jl")

  ## Utilities
  include("utilities/pathways.jl")
  include("utilities/states.jl")
  include("utilities/plotting.jl")

  ## Simulation
  include("simulation/initialization.jl")
  include("simulation/run.jl")
  include("simulation/observe.jl")

  ## Inference
  include("inference/priors.jl")
  include("inference/augmentation.jl")
  include("inference/likelihoods.jl")
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
