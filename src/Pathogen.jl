module Pathogen

  using DataFrames
  using PhyloTrees
  using Distributions
  using RecipesBase

  # New methods
  import
    Base.push!,
    Base.append!,
    Base.rand

  include("core/risks.jl")
  include("core/rates.jl")
  include("core/events.jl")

  include("utilities/pathways.jl")
  include("utilities/states.jl")

  include("simulation/simulate.jl")

  include("inference/priors.jl")
  include("inference/augmentation.jl")
  include("inference/likelihoods.jl")
  include("inference/mcmc.jl")

  include("plotting/helpers.jl")
  include("plotting/recipes.jl")

  # New types and functions
  export
    RiskFunctions,
    RiskParameters,
    Rates,
    Events,
    pathwayto,
    pathwaysto,
    pathwayfrom,
    pathwaysfrom,
    update_events!,
    update_rates!,
    initialize_rates,
    initialize_simulation,
    simulate!,
    generatetree,
    findstate,
    RiskPriors,
    PathogenTrace,
    PathogenIteration,
    PathogenProposal

end # module
