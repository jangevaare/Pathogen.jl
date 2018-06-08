__precompile__()

module Pathogen

  # Dependencies
  using DataFrames
  using Distributions
  # using RecipesBase
  # using ProgressMeter

  import Base.length, Base.convert, Base.show, Base.copy

  # Types
  include("types/EpidemicModel.jl")
  include("types/DiseaseState.jl")
  include("types/Events/Event.jl")
  include("types/Events/EventExtents.jl")
  include("types/Events/EventObservations.jl")
  include("types/Events/EventRates.jl")
  include("types/Risks/RiskFunctions.jl")
  include("types/Risks/RiskParameters.jl")
  include("types/Risks/RiskPriors.jl")
  include("types/Transmissions/Transmission.jl")
  include("types/Transmissions/TransmissionNetwork.jl")
  include("types/Transmissions/TransmissionRates.jl")
  include("types/Simulation.jl")
  include("types/Trace.jl")
  include("types/MCMC.jl")

  # Functions
  include("functions/initialize.jl")
  include("functions/generate.jl")
  include("functions/observe.jl")
  include("functions/update!.jl")
  include("functions/next!.jl")
  include("functions/simulate!.jl")
  include("functions/loglikelihood.jl")

  # Helpers
  include("helpers/RiskFunctions.jl")

  export
    SEIR, SEI, SIR, SI,
    DiseaseState, DiseaseStates,
    RiskFunctions, RiskParameters, RiskPriors,
    Simulation,
    next!, simulate!,
    EventObservations, EventExtents,
    observe,
    MCMC

end
