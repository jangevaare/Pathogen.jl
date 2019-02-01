function initialize(::Type{TransmissionRates},
                    states::Vector{DiseaseState},
                    pop::Population,
                    rf::RiskFunctions{T},
                    rp::RiskParameters{T}) where T <: EpidemicModel
  n_ids = length(states)
  tr = TransmissionRates(n_ids)
  for i in findall(states .== Ref(State_S))
    # External exposure
    #tr.external[i] = rf.susceptibility(rp.susceptibility, pop, i) * rf.sparks(rp.sparks, pop, i)
    tr.external[i] = rf.sparks(rp.sparks, pop, i)
    # Internal exposure
    for k in findall(states .== Ref(State_I))
      tr.internal[k, i] = rf.susceptibility(rp.susceptibility, pop, i) *
                          rf.infectivity(rp.infectivity, pop, i, k) *
                          rf.transmissibility(rp.transmissibility, pop, k)
    end
  end
  @debug "Initialization of $T TransmissionRates complete" external = tr.external ∑external = sum(tr.external) internal = tr.internal ∑internal = sum(tr.internal)
  return tr
end

function initialize(::Type{EventRates},
                    tr::TransmissionRates,
                    states::Vector{DiseaseState},
                    pop::Population,
                    rf::RiskFunctions{T},
                    rp::RiskParameters{T}) where T <: EpidemicModel
  n_ids = length(states)
  rates = EventRates{T}(n_ids)
  for i = 1:n_ids
    if states[i] == State_S
      if T in [SEIR; SEI]
        rates.exposure[i] = tr.external[i] + sum(tr.internal[:,i])
      elseif T in [SIR; SI]
        rates.infection[i] = tr.external[i] + sum(tr.internal[:,i])
      end
    elseif states[i] == State_E
      rates.infection[i] = rf.latency(rp.latency, pop, i)
    elseif states[i] == State_I
      if T in [SEIR; SIR]
        rates.removal[i] = rf.removal(rp.removal, pop, i)
      end
    end
  end
  @debug "Initialization of $T EventRates complete" rates = rates[_state_progressions[T][2:end]]
  return rates
end

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
    events = generate(Events, mcmc)
    rparams = generate(RiskParameters, mcmc)
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
    events = generate(Events, mcmc)
    rparams = generate(RiskParameters, mcmc)
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
