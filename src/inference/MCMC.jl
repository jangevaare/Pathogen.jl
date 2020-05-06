mutable struct MCMC{T <: DiseaseStateSequence}
  event_observations::EventObservations{T}
  event_extents::EventExtents{T}
  population::Population
  starting_states::DiseaseStates
  risk_functions::RiskFunctions{T}
  risk_priors::RiskPriors{T}
  transmission_network_prior::Union{Nothing, TNPrior}
  markov_chains::Vector{MarkovChain{T}}

  function MCMC(obs::EventObservations{T},
                ee::EventExtents{T},
                pop::Population,
                states::DiseaseStates,
                rf::RiskFunctions{T},
                rp::RiskPriors{T};
                tnprior::Union{Nothing, TNPrior}=nothing) where T <: DiseaseStateSequence
    return new{T}(obs, ee, pop, states, rf, rp, tnprior, MarkovChain{T}[])
  end
end

function MCMC(obs::EventObservations{T},
              ee::EventExtents{T},
              pop::Population,
              rf::RiskFunctions{T},
              rp::RiskPriors{T};
              tnprior::Union{Nothing, TNPrior}=nothing) where T <: DiseaseStateSequence
  return MCMC(obs, ee, pop, fill(State_S, individuals(pop)), rf, rp, tnprior=tnprior)
end

function Base.show(io::IO, x::MCMC{T}) where T <: DiseaseStateSequence
  return print(io, "$T model MCMC with $(length(x.markov_chains)) chains")
end

function TransmissionNetworkDistribution(x::MCMC; burnin::Int64=0, thin::Int64=1)
  return TNDistribution([TNDistribution(y, burnin=burnin, thin=thin) for y in x.markov_chains])
end