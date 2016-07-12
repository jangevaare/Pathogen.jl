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

  include("risks.jl")
  include("rates.jl")
  include("events.jl")
  include("pathways.jl")
  include("simulate.jl")
  include("utilities.jl")
  include("plothelpers.jl")
  include("plotrecipes.jl")
  include("infer.jl")

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
