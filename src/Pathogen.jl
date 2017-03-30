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
  include("types/RiskParameterPriors.jl")
  include("types/States.jl")
  include("types/Rates.jl")
  include("types/Network.jl")
  include("types/Events.jl")
  include("types/Event.jl")
  include("types/EventObservations.jl")
  include("types/PathogenIteration.jl")
  include("types/PathogenTrace.jl")

  # Functions
  include("functions/initialize_rates.jl")
  include("functions/generate_event.jl")
  include("functions/update_states!.jl")
  include("functions/update_rates!.jl")
  include("functions/update_events!.jl")
  include("functions/update_network!.jl")
  include("functions/initialize_simulation.jl")
  include("functions/simulate!.jl")
  include("functions/observe.jl")
  include("functions/generate_events.jl")
  include("functions/pathwayto.jl")
  include("functions/pathwayfrom.jl")
  include("functions/generate_tree.jl")
  include("functions/findstate.jl")
  include("functions/popplot.jl")
  include("functions/pathplot.jl")
  include("functions/epiplot.jl")
  include("functions/plot.jl")
  include("functions/propose.jl")
  include("functions/logprior.jl")
  include("functions/loglikelihood.jl")
  include("functions/MHaccept.jl")
  include("initialize_mcmc.jl")
  include("mcmc!.jl")

  export
    SEIR_RiskFunctions,
    SIR_RiskFunctions,
    SEI_RiskFunctions,
    SI_RiskFunctions,
    SEIR_RiskParameters,
    SIR_RiskParameters,
    SEI_RiskParameters,
    SI_RiskParameters,
    SEIR_RiskParameterPriors,
    SIR_RiskParameterPriors,
    SEI_RiskParameterPriors,
    SI_RiskParameterPriors,
    SEIR_EventObservations,
    SIR_EventObservations,
    SEI_EventObservations,
    SI_EventObservations,
    initialize_simulation,
    simulate!,
    observe,
    generate_tree,
    generate_events,
    initialize_mcmc,
    mcmc!

end
