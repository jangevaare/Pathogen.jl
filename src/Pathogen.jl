module Pathogen

  # Dependencies
  using Distributed,
        DataFrames,
        Distributions,
        RecipesBase,
        Logging,
        StatsBase,
        Statistics,
        ProgressMeter,
        LinearAlgebra,
        OnlineStats

  import StatsBase.mean,
         ProgressMeter.next!

  # Types
  include("types/EpidemicModel.jl")
  include("types/DiseaseState.jl")
  include("types/Population.jl")
  include("types/Events/Event.jl")
  include("types/Events/Events.jl")
  include("types/Events/EventExtents.jl")
  include("types/Events/EventObservations.jl")
  include("types/Events/EventRates.jl")
  include("types/Risks/AbstractRisk.jl")
  include("types/Risks/RiskFunctions.jl")
  include("types/Risks/RiskParameters.jl")
  include("types/Risks/RiskPriors.jl")
  include("types/Transmissions/Transmission.jl")
  include("types/Transmissions/TransmissionNetwork.jl")
  include("types/Transmissions/TransmissionRates.jl")
  include("types/Simulation.jl")
  include("types/MarkovChain.jl")
  include("types/MCMC.jl")
  include("types/Transmissions/TransmissionNetworkDistribution.jl")

  # Functions
  include("functions/_pathway_from.jl")
  include("functions/_pathway_to.jl")
  include("functions/_accept.jl")
  include("functions/_bounds.jl")
  include("functions/generate.jl")
  include("functions/initialize.jl")
  include("functions/observe.jl")
  include("functions/update!.jl")
  include("functions/simulate!.jl")
  include("functions/loglikelihood.jl")
  include("functions/logpriors.jl")
  include("functions/start!.jl")
  include("functions/iterate!.jl")

  # Visualization
  include("visualization/epidemic_curve.jl")
  include("visualization/population.jl")
  include("visualization/epidemic_curve.jl")
  include("visualization/transmission_network.jl")
  include("visualization/transmission_network_distribution.jl")
  include("visualization/trace.jl")

  export
    EpidemicModel,
    SEIR, SEI, SIR, SI,
    DiseaseState,
    State_S, State_E, State_I, State_R,
    RiskFunctions, RiskParameters, RiskPriors,
    Population,
    TransmissionNetwork, 
    TransmissionNetworkDistribution, TNDistribution,
    TransmissionNetworkPrior, TNPrior,
    TransmissionNetworkPosterior, TNPosterior,
    Simulation,
    next!, simulate!,
    Events, EventObservations, EventExtents,
    observe,
    MCMC, start!, iterate!
end
