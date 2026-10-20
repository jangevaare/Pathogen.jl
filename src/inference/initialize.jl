function _initialization_attempt(rng::AbstractRNG, mcmc::MCMC{S, M},
                                 max_lposterior::Float64) where {
                                 S <: DiseaseStateSequence,
                                 M <: ILM}
  events = generate(rng, Events, mcmc)
  rparams = generate(rng, RiskParameters, mcmc)
  lprior = logprior(rparams, mcmc.risk_priors)
  if M <: PhyloILM
    sm = generate(rng, mcmc.substitution_model, mcmc.substitution_model_prior)
    lprior += logprior(sm, mcmc.substitution_model_prior)
  end
  llikelihood, tr, tx = loglikelihood(rparams, mcmc.risk_functions, events,
                                      mcmc.population, mcmc.starting_states)
  base_posterior = llikelihood + lprior
  # Network and sequence log probabilities cannot increase this upper bound.
  base_posterior > max_lposterior || return nothing
  _has_network_support(tr, mcmc.transmission_network_prior, tx) || return nothing
  network = generate(rng, TransmissionNetwork, tr, mcmc, tx)
  terms = _network_logprobabilities(network, tr, mcmc.transmission_network_prior,
                                    tx, mcmc.starting_states)
  base_posterior + sum(terms) > max_lposterior || return nothing
  if M <: PhyloILM
    tree, obs_nodes = generate(Tree, events, mcmc.event_observations, network)
    leaf_data = Dict{Int64, RNASeq}()
    for k in eachindex(obs_nodes)
      if !isnothing(obs_nodes[k])
        leaf_data[obs_nodes[k]] = mcmc.event_observations.seq[k]
      end
    end
    base_posterior += PhyloModels.loglikelihood(tree, sm, leaf_data)
  end
  base_posterior + sum(terms) > max_lposterior || return nothing
  if M <: PhyloILM
    return MarkovChain{S, M}(events, network, rparams, sm, base_posterior;
                             network_logprobabilities=terms)
  else
    return MarkovChain{S, M}(events, network, rparams, base_posterior;
                             network_logprobabilities=terms)
  end
end

function initialize(rng::AbstractRNG, ::Type{MarkovChain},
                    mcmc::MCMC{S, M},
                    progress_channel::RemoteChannel;
                    attempts::Int64=1000) where {
                    S <: DiseaseStateSequence,
                    M <: ILM}
  if attempts <= 0
    throw(ErrorException("Must have at least 1 initialization attempt"))
  end
  local markov_chain
  max_lposterior = -Inf
  for i in 1:attempts
    initialization = _initialization_attempt(rng, mcmc, max_lposterior)
    if initialization !== nothing
      markov_chain = initialization
      max_lposterior = markov_chain.log_posterior[1]
    end
    put!(progress_channel, true)
  end
  if max_lposterior == -Inf
    @error "Failed to initialize Markov Chain"
    markov_chain = nothing
  end
  return markov_chain
end


function initialize(rng::AbstractRNG, ::Type{MarkovChain},
                    mcmc::MCMC{S, M};
                    attempts::Int64=1000) where {
                    S <: DiseaseStateSequence,
                    M <: ILM}
  if attempts <= 0
    throw(ErrorException("Must have at least 1 initialization attempt"))
  end
  max_lposterior = -Inf
  local markov_chain
  pmeter = Progress(attempts, "Initialization progress")
  for i in 1:attempts
    initialization = _initialization_attempt(rng, mcmc, max_lposterior)
    if initialization !== nothing
      markov_chain = initialization
      max_lposterior = markov_chain.log_posterior[1]
    end
    next!(pmeter)
  end
  finish!(pmeter)
  if max_lposterior == -Inf
    @error "Failed to initialize Markov Chain"
    markov_chain = nothing
  end
  return markov_chain
end

initialize(::Type{MarkovChain}, args...; kwargs...) = initialize(_default_rng(), MarkovChain, args...; kwargs...)
_initialization_attempt(mcmc::MCMC, max_lposterior::Float64) = _initialization_attempt(_default_rng(), mcmc, max_lposterior)
