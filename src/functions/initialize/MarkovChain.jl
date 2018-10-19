function initialize(::Type{MarkovChain},
                    mcmc::MCMC{T},
                    progress_channel::RemoteChannel;
                    attempts::Int64=1000) where T <: EpidemicModel
  if attempts <= 0
    @error "Must have at least 1 initialization attempt"
  end
  max_lposterior = -Inf
  local markov_chain
  for i in 1:attempts
    rparams = generate(RiskParameters, mcmc)
    events = generate(Events, rparams, mcmc)
    lprior = logpriors(rparams, mcmc.risk_priors)
    llikelihood, network = loglikelihood(rparams,
                                         mcmc.risk_functions,
                                         events,
                                         mcmc.population,
                                         mcmc.starting_states,
                                         early_decision_value = max_lposterior - lprior)
    lposterior = llikelihood + lprior
    if lposterior > max_lposterior
      markov_chain = MarkovChain(events, network, rparams, lposterior)
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
                    mcmc::MCMC{T};
                    attempts::Int64=1000) where T <: EpidemicModel
  if attempts <= 0
    @error "Must have at least 1 initialization attempt"
  end
  max_lposterior = -Inf
  local markov_chain
  pmeter = Progress(attempts, "Initialization progress")
  for i in 1:attempts
    @debug "Beginning MarkovChain initialization attempt $i"
    next!(pmeter)
    rparams = generate(RiskParameters, mcmc)
    events = generate(Events, rparams, mcmc)
    lprior = logpriors(rparams, mcmc.risk_priors)
    llikelihood, network = loglikelihood(rparams,
                                         mcmc.risk_functions,
                                         events,
                                         mcmc.population,
                                         mcmc.starting_states,
                                         early_decision_value = max_lposterior - lprior)
    lposterior = llikelihood + lprior
    if lposterior > max_lposterior
      markov_chain = MarkovChain(events, network, rparams, lposterior)
      max_lposterior = lposterior
    end
  end
  finish!(pmeter)
  if max_lposterior == -Inf
    @error "Failed to initialize Markov Chain"
  end
  return markov_chain
end
