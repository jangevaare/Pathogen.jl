"""
Initiate an Trace object by sampling from specified prior distributions
"""
function initialize(ilm_priors::SIR_priors,
                    mutation_priors::JC69_priors,
                    detection_priors::Lag_priors,
                    obs::SIR_observed,
                    limit=500::Int,
                    debug=false::Bool,
                    dist=Euclidean())
  count = 0
  lp = -Inf

  # Retry initialization until non-negative infinity loglikelihood
  while (lp == -Inf || lp == Inf) && count < limit
    count += 1
    ilm_params = randprior(ilm_priors)
    mutation_params = randprior(mutation_priors)
    detection_params = randprior(detection_priors)
    aug = propose_augment(detection_params[1], obs, debug)
    lp, network_rates = SIR_loglikelihood(ilm_params[1],
                                          ilm_params[2],
                                          ilm_params[3],
                                          ilm_params[4],
                                          aug,
                                          obs,
                                          debug,
                                          dist)
    lp += logprior(ilm_priors,
                   ilm_params,
                   debug)

    lp += logprior(detection_priors,
                   detection_params,
                   debug)
    if Inf > lp > -Inf
      network = propose_network(network_rates,
                                debug)
      lp += exposure_network_loglikelihood(network,
                                           network_rates,
                                           debug)
      lp += detection_loglikelihood(detection_params[1],
                                    aug,
                                    network,
                                    obs,
                                    debug)
    end
    if Inf > lp > -Inf
      lp += logprior(mutation_priors,
                     mutation_params,
                     debug)
    end
    if Inf > lp > -Inf
      lp += phylogenetic_network_loglikelihood(obs,
                                               aug,
                                               network,
                                               jc69p([mutation_params[1]]),
                                               debug)
    end
    if Inf > lp > -Inf
      return SIR_trace([ilm_params[1]],
                       [ilm_params[2]],
                       [ilm_params[3]],
                       [ilm_params[4]],
                       [aug],
                       Array[network_rates],
                       Array[network],
                       [lp]),
              Lag_trace([detection_params[1]]),
              JC69_trace([mutation_params[1]])
    end
    @assert(count < limit && Inf > lp > -Inf, "Failed to initialize")
  end
end


"""
Initiate an Trace object by sampling from specified prior distributions
"""
function initialize(ilm_priors::SIR_priors,
                    detection_priors::Lag_priors,
                    obs::SEIR_observed,
                    limit=500::Int,
                    debug=false::Bool,
                    dist=Euclidean())
  count = 0
  lp = -Inf

  # Retry initialization until non-negative infinity loglikelihood
  while (lp == -Inf || lp == Inf) && count < limit
    count += 1
    ilm_params = randprior(ilm_priors)
    detection_params = randprior(detection_priors)
    aug = propose_augment(detection_params[1],
                          obs,
                          debug)
    lp, network_rates = SIR_loglikelihood(ilm_params[1],
                                          ilm_params[2],
                                          ilm_params[3],
                                          ilm_params[4],
                                          aug,
                                          obs,
                                          debug,
                                          dist)
    lp += logprior(ilm_priors, ilm_params, debug)
    lp += logprior(detection_priors, detection_params)
    if Inf > lp > -Inf
      network = propose_network(network_rates, debug)
      lp += detection_loglikelihood(detection_params[1],
                                   aug,
                                   network,
                                   obs,
                                   debug)
      lp += exposure_network_loglikelihood(network, network_rates, debug)
    end
    if Inf > lp > -Inf
      return SIR_trace([ilm_params[1]],
                       [ilm_params[2]],
                       [ilm_params[3]],
                       [ilm_params[4]],
                       [aug],
                       Array[network_rates],
                       Array[network],
                       [lp]),
              Lag_trace([detection_params[1]])
    end
    @assert(count < limit && Inf > lp > -Inf, "Failed to initialize")
  end
end


"""
Initiate an Trace object by sampling from specified prior distributions
"""
function initialize(ilm_priors::SIR_priors,
                    obs::SIR_observed,
                    limit=500::Int,
                    debug=false::Bool,
                    dist=Euclidean())
  count = 0
  lp = -Inf

  # Retry initialization until non-negative infinity loglikelihood
  while (lp == -Inf || lp == Inf) && count < limit
    count += 1
    ilm_params = randprior(ilm_priors)
    aug = propose_augment(obs,
                          debug)
    lp, network_rates = SIR_loglikelihood(ilm_params[1],
                                          ilm_params[2],
                                          ilm_params[3],
                                          ilm_params[4],
                                          aug,
                                          obs,
                                          debug,
                                          dist)
    lp += logprior(ilm_priors, ilm_params, debug)
    if Inf > lp > -Inf
      network = propose_network(network_rates, debug)
      lp += exposure_network_loglikelihood(network, network_rates, debug)
    end
    if Inf > lp > -Inf
      return SIR_trace([ilm_params[1]],
                       [ilm_params[2]],
                       [ilm_params[3]],
                       [ilm_params[4]],
                       [aug],
                       Array[network_rates],
                       Array[network],
                       [lp]),
    end
    @assert(count < limit && Inf > lp > -Inf, "Failed to initialize")
  end
end
