function initialize(::Type{MarkovChain},
                    mcmc::MCMC{T},
                    progress_channel::RemoteChannel;
                    attempts::Int64=1000) where T <: EpidemicModel
  if attempts <= 0
    error("Must have at least 1 initialization attempt")
  end
  bested = 0
  max_lposterior = -Inf
  local markov_chain
  for i in 1:attempts
    events = generate(Events, mcmc)
    rparams = generate(RiskParameters, mcmc)
    @debug "Risk function parameters for attempt $i" θ = convert(Vector{Float64}, rparams)
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
      markov_chain = MarkovChain(events, network, rparams, lposterior)
      max_lposterior = lposterior
      bested += 1
    end
    put!(progress_channel, true)
  end
  if max_lposterior == -Inf
    error("Failed to initialize Markov Chain")
  end
  @debug "Initialization improved upon $bested times"
  return markov_chain
end

function initialize(::Type{MarkovChain},
                    mcmc::MCMC{T};
                    attempts::Int64=1000) where T <: EpidemicModel
  if attempts <= 0
    error("Must have at least 1 initialization attempt")
  end
  max_lposterior = -Inf
  local markov_chain
  bested = 0
  pmeter = Progress(attempts, "Initialization progress")
  for i in 1:attempts
    next!(pmeter)
    events = generate(Events, mcmc)
    rparams = generate(RiskParameters, mcmc)
    @debug "Risk function parameters for attempt $i" θ = convert(Vector{Float64}, rparams)
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
      markov_chain = MarkovChain(events, network, rparams, lposterior)
      max_lposterior = lposterior
      bested += 1
    end
  end
  finish!(pmeter)
  if max_lposterior == -Inf
    error("Failed to initialize Markov Chain")
    markov_chain = nothing
  end
  @debug "Initialization improved upon $bested times"
  return markov_chain
end
