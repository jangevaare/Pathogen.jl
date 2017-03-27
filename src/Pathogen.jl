__precompile__()

module Pathogen

  # Dependencies
  using DataFrames
  using PhyloTrees
  using PhyloModels
  using Distributions
  using RecipesBase
  using ProgressMeter

  # Methods
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

  # Types
  include("types/RiskFunctions.jl")
  include("types/RiskParameters.jl")
  include("types/States.jl")
  include("types/Rates.jl")
  include("types/Network.jl")
  include("types/Events.jl")
  include("types/Event.jl")
  include("types/EventObservations.jl")

  # Functions
  include("functions/initialize_rates.jl")
  include("functions/generate_event.jl")
  include("functions/update_states!.jl")
  include("functions/update_rates!.jl")
  include("functions/update_events!.jl")
  include("functions/update_network!.jl")
  include("functions/initialize_simulation.jl")
  include("functions/simulate!.jl")
  include("functions/findstate.jl")

  export
    SEIR_RiskFunctions,
    SIR_RiskFunctions,
    SEI_RiskFunctions,
    SI_RiskFunctions,
    SEIR_RiskParameters,
    SIR_RiskParameters,
    SEI_RiskParameters,
    SI_RiskParameters,
    initialize_simulation,
    simulate!

end
