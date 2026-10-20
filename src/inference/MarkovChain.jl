mutable struct MarkovChain{S <: DiseaseStateSequence, M <: ILM}
  iterations::Int64
  events::Vector{Events{S}}
  transmission_network::Vector{TransmissionNetwork}
  risk_parameters::Vector{RiskParameters{S}}
  substitution_model::Union{Nothing, Vector{NucleicAcidSubstitutionModel}}
  log_posterior::Vector{Float64}
  # The base includes epidemic likelihood, parameter priors and any sequence
  # likelihood. Retain parent probabilities separately for proposal cancellation.
  log_base_posterior::Float64
  network_logprobabilities::Vector{Float64}
  Σrp::OnlineStats.CovMatrix
  Σsm::Union{Nothing, OnlineStats.CovMatrix}

  function MarkovChain{S, M}(e::Events{S},
                             n::TransmissionNetwork,
                             rp::RiskParameters{S},
                             lp::Float64;
                             network_logprobabilities::Vector{Float64}=Float64[]) where {
                             S <: DiseaseStateSequence,
                             M <: TNILM}
    return new{S, M}(0, [e], [n], [rp], nothing, [lp + sum(network_logprobabilities)],
                     lp, network_logprobabilities, OnlineStats.CovMatrix(), nothing)
  end

  function MarkovChain{S, M}(e::Events{S},
                             n::TransmissionNetwork,
                             rp::RiskParameters{S},
                             sm::NucleicAcidSubstitutionModel,
                             lp::Float64;
                             network_logprobabilities::Vector{Float64}=Float64[]) where {
                             S <: DiseaseStateSequence,
                             M <: PhyloILM}
    return new{S, M}(0, [e], [n], [rp], [sm], [lp + sum(network_logprobabilities)],
                     lp, network_logprobabilities, OnlineStats.CovMatrix(), OnlineStats.CovMatrix())
  end
end

function Base.length(x::MarkovChain)
  return x.iterations
end

function Base.show(io::IO, x::MarkovChain{S, M}) where {
                   S <: DiseaseStateSequence,
                   M <: ILM}
  return print(io, "$S $M Markov chain (iterations = $(length(x)))")
end

function TransmissionNetworkDistribution(x::MarkovChain; burnin::Int64=0, thin::Int64=1)
  return TNDistribution(x.transmission_network[1+burnin:thin:end])
end

# Direct constructors accept the historical score without parent probabilities.
# If initialization did not supply them, compute them on first use.
function _ensure_network_probabilities!(mc::MarkovChain, mcmc, work, cache)
  if isempty(mc.network_logprobabilities)
    _, rates, transmissions = loglikelihood(mc.risk_parameters[end], mcmc.risk_functions,
      mc.events[end], mcmc.population, mcmc.starting_states;
      loglikelihood_output=false, transmission_rate_cache=cache, workspace=work)
    terms = _network_logprobabilities(mc.transmission_network[end], rates,
      mcmc.transmission_network_prior, transmissions, mcmc.starting_states)
    isfinite(sum(terms)) || throw(ArgumentError("The starting network has zero conditional probability"))
    mc.network_logprobabilities = terms
    mc.log_posterior[end] = mc.log_base_posterior + sum(terms)
  end
  return mc.network_logprobabilities
end
