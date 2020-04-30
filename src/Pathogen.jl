module Pathogen

  # Dependencies
  using Distributed,
        DataFrames,
        Distributions,
        RecipesBase,
        Logging,
        StaticArrays,
        StatsBase,
        Statistics,
        ProgressMeter,
        LinearAlgebra,
        OnlineStats,
        PhyloModels

  import StatsBase.mean,
         ProgressMeter.next!,
         Base.summary

  # Core
  include("core/IndividualLevelModel.jl")
  include("core/DiseaseStateSequence.jl")
  include("core/DiseaseState.jl")
  include("core/Population.jl")
  include("core/Events/Event.jl")
  include("core/Events/Events.jl")
  include("core/Events/EventRates.jl")
  include("core/Events/EventObservations.jl")
  include("core/Risks/AbstractRisk.jl")
  include("core/Risks/RiskFunctions.jl")
  include("core/Risks/RiskParameters.jl")
  include("core/Transmissions/Transmission.jl")
  include("core/Transmissions/AbstractTransmissionNetwork.jl")
  include("core/Transmissions/TransmissionNetwork.jl")
  include("core/Transmissions/TransmissionRates.jl")
  include("core/Transmissions/TransmissionRateCache.jl")
  include("core/initialize/TransmissionRates.jl")
  include("core/initialize/EventRates.jl")
  include("core/update!/DiseaseStates.jl")
  include("core/update!/Events.jl")
  include("core/update!/TransmissionNetwork.jl")
  include("core/update!/EventRates.jl")
  include("core/update!/TransmissionRates.jl")

  # Simulation
  include("simulation/Simulation.jl")
  include("simulation/generate/Event.jl")
  include("simulation/generate/Transmission.jl")
  include("simulation/generate/Tree.jl")
  include("simulation/simulate!.jl")
  include("simulation/update!.jl")
  include("simulation/generate/EventObservations.jl")

  # Inference
  include("inference/TransmissionNetworkDistribution.jl")
  include("inference/RiskPriors.jl")
  include("inference/EventExtents.jl")
  include("inference/MarkovChain.jl")
  include("inference/MCMC.jl")
  include("inference/generate/Event.jl")
  include("inference/generate/Events.jl")
  include("inference/generate/NucleicAcidSubstitutionModel.jl")
  include("inference/generate/RiskParameters.jl")
  include("inference/generate/Transmission.jl")
  include("inference/generate/TransmissionNetwork.jl")
  include("inference/generate/Tree.jl")
  include("inference/loglikelihood.jl")
  include("inference/initialize.jl")
  include("inference/update!.jl")
  include("inference/logprior.jl")
  include("inference/start!.jl")
  include("inference/iterate!.jl")
  include("inference/summary.jl")

  # Visualization
  include("visualization/epidemic_curve.jl")
  include("visualization/observations.jl")
  include("visualization/population.jl")
  include("visualization/transmission_network.jl")
  include("visualization/transmission_network_distribution.jl")
  include("visualization/trace.jl")
  include("visualization/out_degree.jl")

  # re-export PhyloModels, DataFrames, Distributions
  for reexport_pkg in [PhyloModels, DataFrames, Distributions]
    for name in names(reexport_pkg)
      @eval export $(name)
    end
  end

  export
    IndividualLevelModel, ILM,
    PhylodynamicILM, PhyloILM,
    TransmissionNetworkILM, TNILM,
    DiseaseStateSequence,
    SEIR, SEI, SIR, SI,
    DiseaseState, DiseaseStates,
    State_S, State_E, State_I, State_R,
    RiskFunctions, RiskParameters, RiskPriors,
    Population,
    TransmissionNetwork,
    TransmissionNetworkDistribution, TNDistribution,
    TransmissionNetworkPrior, TNPrior,
    TransmissionNetworkPosterior, TNPosterior,
    Transmission,
    EndogenousTransmission, ExogenousTransmission, NoTransmission,
    individuals,
    Simulation, simulate!, generate,
    Events, EventObservations, EventExtents,
    MCMC, start!, iterate!, summary
end
