mutable struct MCMC{T <: EpidemicModel}
  event_observations::EventObservations{T}
  event_extents::EventExtents{T}
  population::Population
  starting_states::Vector{DiseaseState}
  risk_functions::RiskFunctions{T}
  risk_priors::RiskPriors{T}
  transmission_network_prior::Union{Nothing, TNPrior}
  markov_chains::Vector{MarkovChain{T}}

  function MCMC(obs::EventObservations{T},
                ee::EventExtents{T},
                pop::Population,
                states::Vector{DiseaseState},
                rf::RiskFunctions{T},
                rp::RiskPriors{T};
                tnprior::Union{Nothing, TNPrior}=nothing) where T <: EpidemicModel
    return new{T}(obs, ee, pop, states, rf, rp, tnprior, MarkovChain{T}[])
  end
end

function MCMC(obs::EventObservations{T},
              ee::EventExtents{T},
              pop::Population,
              rf::RiskFunctions{T},
              rp::RiskPriors{T};
              tnprior::Union{Nothing, TNPrior}=nothing) where T <: EpidemicModel
  return MCMC(obs, ee, pop, fill(State_S, pop.individuals), rf, rp, tnprior=tnprior)
end

function Base.show(io::IO, x::MCMC{T}) where T <: EpidemicModel
  return print(io, "$T model MCMC with $(length(x.markov_chains)) chains")
end

function TransmissionNetworkDistribution(x::MCMC)
  return TNDistribution([TNDistribution(y) for y in x.markov_chains])
end