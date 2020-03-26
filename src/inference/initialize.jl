function initialize(::Type{MarkovChain},
                    mcmc::MCMC{S, M},
                    progress_channel::RemoteChannel;
                    attempts::Int64=1000) where {
                    S <: DiseaseStateSequence,
                    M <: TNILM}
  if attempts <= 0
    @error "Must have at least 1 initialization attempt"
  end
  max_lposterior = -Inf
  local markov_chain
  for i in 1:attempts
    events = generate(Events, mcmc)
    rparams = generate(RiskParameters, mcmc)
    lprior = logpriors(rparams, mcmc.risk_priors)
    llikelihood, tr, tx = loglikelihood(rparams,
                                        mcmc.risk_functions,
                                        events,
                                        mcmc.population,
                                        mcmc.starting_states,
                                        early_decision_value = max_lposterior - lprior)
    lposterior = llikelihood + lprior
    if lposterior > max_lposterior
      network = generate(TransmissionNetwork, tr, mcmc, tx)
      markov_chain = MarkovChain{M}(events, network, rparams, lposterior)
      max_lposterior = lposterior
    end
    put!(progress_channel, true)
  end
  if max_lposterior == -Inf
    @error "Failed to initialize Markov Chain"
  end
  return markov_chain
end

function initialize(::Type{MarkovChain},
                    mcmc::MCMC{S, M};
                    attempts::Int64=1000) where {
                    S <: DiseaseStateSequence,
                    M <: TNILM}
  if attempts <= 0
    @error "Must have at least 1 initialization attempt"
  end
  max_lposterior = -Inf
  local markov_chain
  pmeter = Progress(attempts, "Initialization progress")
  for i in 1:attempts
    @debug "Beginning MarkovChain initialization attempt $i"
    next!(pmeter)
    events = generate(Events, mcmc)
    rparams = generate(RiskParameters, mcmc)
    lprior = logpriors(rparams, mcmc.risk_priors)
    llikelihood, tr, tx = loglikelihood(rparams,
                                        mcmc.risk_functions,
                                        events,
                                        mcmc.population,
                                        mcmc.starting_states,
                                        early_decision_value = max_lposterior - lprior)
    lposterior = llikelihood + lprior
    if lposterior > max_lposterior
      network = generate(TransmissionNetwork, tr, mcmc, tx)
      markov_chain = MarkovChain{M}(events, network, rparams, lposterior)
      max_lposterior = lposterior
    end
  end
  finish!(pmeter)
  if max_lposterior == -Inf
    @error "Failed to initialize Markov Chain"
    markov_chain = nothing
  end
  return markov_chain
end
