mutable struct MCMC{S <: DiseaseStateSequence, M <: ILM}
  event_observations::EventObservations{S, M}
  event_extents::EventExtents{S}
  population::Population
  starting_states::DiseaseStates
  risk_functions::RiskFunctions{S}
  risk_priors::RiskPriors{S}
  substitution_model::Union{Nothing, Type{<:NucleicAcidSubstitutionModel}}
  substitution_model_priors::Union{Nothing, Vector{UnivariateDistribution}}
  transmission_network_prior::Union{Nothing, TNPrior}
  markov_chains::Vector{MarkovChain{S}}

  function MCMC{S, M}(obs::EventObservations{S},
                      ee::EventExtents{S},
                      pop::Population,
                      states::DiseaseStates,
                      rf::RiskFunctions{S},
                      rp::RiskPriors{S};
                      tnprior::Union{Nothing, TNPrior}=nothing) where {
                      S <: DiseaseStateSequence,
                      M <: TNILM}
    return new{S, M}(obs, ee, pop, states, rf, rp, nothing, nothing, tnprior, MarkovChain{S, M}[])
  end

  function MCMC{S, M}(obs::EventObservations{S},
                      ee::EventExtents{S},
                      pop::Population,
                      states::DiseaseStates,
                      rf::RiskFunctions{S},
                      rp::RiskPriors{S},
                      sm::Type{N},
                      smp::Vector{UnivariateDistribution};
                      tnprior::Union{Nothing, TNPrior}=nothing) where {
                      S <: DiseaseStateSequence,
                      M <: PhyloILM,
                      N <: NucleicAcidSubstitutionModel}
    return new{S, M}(obs, ee, pop, states, rf, rp, sm, smp, tnprior, MarkovChain{S, M}[])
  end
end

function MCMC(obs::EventObservations{S},
              ee::EventExtents{S},
              pop::Population,
              rf::RiskFunctions{S},
              rp::RiskPriors{S};
              tnprior::Union{Nothing, TNPrior}=nothing) where {
              S <: DiseaseStateSequence,
              M <: TNILM}
  return MCMC{S, M}(obs, ee, pop, fill(State_S, pop.individuals), rf, rp, nothing, nothing, tnprior=tnprior)
end

function MCMC(obs::EventObservations{S},
                 ee::EventExtents{S},
                 pop::Population,
                 rf::RiskFunctions{S},
                 rp::RiskPriors{S},
                 sm::Type{N},
                 smp::Vector{UnivariateDistribution};
                 tnprior::Union{Nothing, TNPrior}=nothing) where {
                 S <: DiseaseStateSequence,
                 M <: PhyloILM,
                 N <: NucleicAcidSubstitutionModel}
  return MCMC{S, M}(obs, ee, pop, fill(State_S, pop.individuals), rf, rp, sm, smp, tnprior, MarkovChain{S}[])
end

function Base.show(io::IO, x::MCMC{S, M}) where {
                   S <: DiseaseStateSequence,
                   M <: ILM}
  return print(io, "$S $M MCMC with $(length(x.markov_chains)) Markov chains")
end

function TransmissionNetworkDistribution(x::MCMC)
  return TNDistribution([TNDistribution(y) for y in x.markov_chains])
end

function TransmissionNetworkDistribution(iter, x::MCMC)
  return TNDistribution([TNDistribution(iter, y) for y in x.markov_chains])
end