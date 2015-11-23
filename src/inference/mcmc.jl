"""
Calculate the log prior from prior distributions and specified parameter values
"""
function logprior(priors::Priors, params::Vector{Float64}, debug=false::Bool)
  @assert(length(params) == length(fieldnames(priors)),
          "Mismatch between parameter vector and prior")
  lprior = 0.
  for i = 1:length(params)
    lprior += logpdf(priors.(fieldnames(priors)[i]), params[i])
  end
  debug && println("$(typeof(priors)) log prior: $(round(lprior,3))")
  return lprior
end


"""
Randomly generate a parameter vector from specified priors
"""
function randprior(priors::Priors)
  params = [rand(priors.(fieldnames(priors)[1]))]
  for i = 2:length(fieldnames(priors))
    push!(params, rand(priors.(fieldnames(priors)[i])))
  end
  return params
end


"""
Initiate an Trace object by sampling from specified prior distributions
"""
function initialize(ilm_priors::SEIR_priors,
                    mutation_priors::JC69_priors,
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
    mutation_params = randprior(mutation_priors)
    detection_params = randprior(detection_priors)
    aug = propose_augment(ilm_params[4], detection_params[1], obs, debug)
    lp, network_rates = SEIR_loglikelihood(ilm_params[1],
                                           ilm_params[2],
                                           ilm_params[3],
                                           ilm_params[4],
                                           ilm_params[5],
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
      return SEIR_trace([ilm_params[1]],
                        [ilm_params[2]],
                        [ilm_params[3]],
                        [ilm_params[4]],
                        [ilm_params[5]],
                        [aug],
                        Array[network_rates],
                        Array[network], [lp]),
              Lag_trace([detection_params[1]]),
              JC69_trace([mutation_params[1]])
    end
    @assert(count < limit && Inf > lp > -Inf, "Failed to initialize")
  end
end


"""
Initiate an Trace object by sampling from specified prior distributions
"""
function initialize(ilm_priors::SEIR_priors,
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
    aug = propose_augment(ilm_params[4],
                          detection_params[1],
                          obs,
                          debug)
    lp, network_rates = SEIR_loglikelihood(ilm_params[1],
                                           ilm_params[2],
                                           ilm_params[3],
                                           ilm_params[4],
                                           ilm_params[5],
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
      return SEIR_trace([ilm_params[1]],
                        [ilm_params[2]],
                        [ilm_params[3]],
                        [ilm_params[4]],
                        [ilm_params[5]],
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
function initialize(ilm_priors::SEIR_priors,
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
    aug = propose_augment(ilm_params[4],
                          obs,
                          debug)
    lp, network_rates = SEIR_loglikelihood(ilm_params[1],
                                           ilm_params[2],
                                           ilm_params[3],
                                           ilm_params[4],
                                           ilm_params[5],
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
      return SEIR_trace([ilm_params[1]],
                        [ilm_params[2]],
                        [ilm_params[3]],
                        [ilm_params[4]],
                        [ilm_params[5]],
                        [aug],
                        Array[network_rates],
                        Array[network],
                        [lp])
    end
    @assert(count < limit && Inf > lp > -Inf, "Failed to initialize")
  end
end


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


"""
A simple function for Metropolis-Hastings rejection using log posteriors
"""
function MHreject(lp1::Float64, lp2::Float64, debug=false::Bool)
  @assert(lp1 < Inf && lp2 < Inf, "Infinite log posterior detected")
  reject = true
  if lp1 >= lp2
    debug && println("MCMC proposal accepted")
    reject = false
  elseif exp(lp1 - lp2) >= rand()
    debug && println("MCMC proposal probabilistically accepted")
    reject = false
  end
  debug && reject && println("MCMC Proposal rejected")
  return reject
end


function MCMC(n::Int64,
              transition_cov::Array{Float64},
              ilm_trace::SEIR_trace,
              detection_trace::Lag_trace,
              mutation_trace::JC69_trace,
              ilm_priors::SEIR_priors,
              detection_priors::Lag_priors,
              mutation_priors::JC69_priors,
              obs::SEIR_observed,
              debug=false::Bool,
              progress=true::Bool,
              dist=Euclidean())

  @assert(size(transition_cov) == (7,7),
          "Transition kernel covariance matrix must be a positive definite 7x7 matrix")
  progressbar = Progress(n, 5, "Performing $n MCMC iterations...", 25)
  debug && println("MCMC transition kernel covariance matrix:")
  debug && println(round(transition_cov, 3))
  rejects = 0
  for i = 1:n
    progress && !debug && next!(progressbar)
    debug && println("")
    debug && println("Performing the $(i)th MCMC iteration")
    if mod(i, 1) == 0
      step = rand(MvNormal(transition_cov))
    else
      step = fill(0., 7)
    end
    ilm_proposal = [ilm_trace.α[end],
                    ilm_trace.β[end],
                    ilm_trace.η[end],
                    ilm_trace.ρ[end],
                    ilm_trace.γ[end]] .+ step[1:5]
    detection_proposal = [detection_trace.ν[end]] .+ step[6]
    mutation_proposal = [mutation_trace.λ[end]] .+ step[7]

    debug && println("Proposed parameter changes: $(round(step, 3))")
    debug && println("ILM proposal: $(round(ilm_proposal, 3))")
    debug && println("Detection proposal: $(round(detection_proposal, 3))")
    debug && println("Mutation proposal: $(round(mutation_proposal, 3))")
    lp = logprior(ilm_priors, ilm_proposal, debug)
    lp += logprior(detection_priors, detection_proposal, debug)
    lp += logprior(mutation_priors, mutation_proposal, debug)

    if lp > -Inf
      if mod(i, 1) == 0
        # Generate data augmentation proposal
        aug = propose_augment(ilm_proposal[4],
                              detection_proposal[1],
                              ilm_trace.network[end],
                              ilm_trace.aug[end],
                              obs,
                              rand(Poisson(1)),
                              debug)
      else
        aug = ilm_trace.aug[end]
      end
    end
    if lp > -Inf
      # SEIR loglikelihood
      ll, network_rates = SEIR_loglikelihood(ilm_proposal[1],
                                             ilm_proposal[2],
                                             ilm_proposal[3],
                                             ilm_proposal[4],
                                             ilm_proposal[5],
                                             aug,
                                             obs,
                                             debug,
                                             dist)
      lp += ll
    end
    if lp > -Inf
      if mod(i, 1) == 0
        # Generate network proposal
        network = propose_network(network_rates,
                                  ilm_trace.network[end],
                                  debug,
                                  rand(Poisson(1)),
                                  "multinomial")
      else
        network = ilm_trace.network[end]
      end
      lp += exposure_network_loglikelihood(network, network_rates, debug)
    end
    if lp > -Inf
      lp += phylogenetic_network_loglikelihood(obs,
                                               aug,
                                               network,
                                               jc69p(mutation_proposal),
                                               debug)
    end
    if lp > -Inf
      lp += detection_loglikelihood(detection_proposal[1],
                                    aug,
                                    network,
                                    obs,
                                    debug)
    end
    # Acceptance/rejection
    reject = MHreject(lp, ilm_trace.logposterior[end], debug)

    if reject
      aug = ilm_trace.aug[end]
      network_rates = ilm_trace.network_rates[end]
      network = ilm_trace.network[end]
      lp = ilm_trace.logposterior[end]
      ilm_proposal = [ilm_trace.α[end],
                      ilm_trace.β[end],
                      ilm_trace.η[end],
                      ilm_trace.ρ[end],
                      ilm_trace.γ[end]]
      detection_proposal = [detection_trace.ν[end]]
      mutation_proposal = [mutation_trace.λ[end]]
      rejects += 1
    end

    push!(ilm_trace.aug, aug)
    push!(ilm_trace.network_rates, network_rates)
    push!(ilm_trace.network, network)
    push!(ilm_trace.logposterior, lp)
    push!(ilm_trace.α, ilm_proposal[1])
    push!(ilm_trace.β, ilm_proposal[2])
    push!(ilm_trace.η, ilm_proposal[3])
    push!(ilm_trace.ρ, ilm_proposal[4])
    push!(ilm_trace.γ, ilm_proposal[5])
    push!(detection_trace.ν, detection_proposal[1])
    push!(mutation_trace.λ, mutation_proposal[1])
  end
  println("MCMC acceptance rate: $(round(1.0-(rejects/n),4))")
  return ilm_trace, detection_trace, mutation_trace
end


function MCMC(n::Int64,
              transition_cov::Array{Float64},
              ilm_priors::SEIR_priors,
              detection_priors::Lag_priors,
              mutation_priors::JC69_priors,
              obs::SEIR_observed,
              debug=false::Bool,
              progress=true::Bool,
              dist=Euclidean(),
              init_limit=500::Int64)

  ilm_trace, detection_trace, mutation_trace = initialize(ilm_priors,
                                                          mutation_priors,
                                                          detection_priors,
                                                          obs,
                                                          init_limit,
                                                          debug,
                                                          dist)

  return MCMC(n,
              transition_cov,
              ilm_trace,
              detection_trace,
              mutation_trace,
              ilm_priors,
              detection_priors,
              mutation_priors,
              obs,
              debug,
              progress,
              dist)

  end


function MCMC(n::Int64,
              ilm_priors::SEIR_priors,
              detection_priors::Lag_priors,
              mutation_priors::JC69_priors,
              obs::SEIR_observed,
              debug=false::Bool,
              progress=true::Bool,
              dist=Euclidean(),
              init_limit=500::Int64)

  transition_cov = diagm([var(ilm_priors.α),
                          var(ilm_priors.β),
                          var(ilm_priors.η),
                          var(ilm_priors.ρ),
                          var(ilm_priors.γ),
                          var(detection_priors.ν),
                          var(mutation_priors.λ)])*(2.38^2)/7.

  return MCMC(n,
              transition_cov,
              ilm_priors,
              detection_priors,
              mutation_priors,
              obs,
              debug,
              progress,
              dist,
              init_limit)
end


function MCMC(n::Int64,
              transition_cov::Array{Float64},
              ilm_trace::SEIR_trace,
              detection_trace::Lag_trace,
              ilm_priors::SEIR_priors,
              detection_priors::Lag_priors,
              obs::SEIR_observed,
              debug=false::Bool,
              progress=true::Bool,
              dist=Euclidean())

  @assert(size(transition_cov) == (6,6),
          "Transition kernel covariance matrix must be a positive definite 6x6 matrix")
  progressbar = Progress(n, 5, "Performing $n MCMC iterations...", 25)
  debug && println("MCMC transition kernel covariance matrix:")
  debug && println(round(transition_cov, 3))
  rejects = 0
  for i = 1:n
    progress && !debug && next!(progressbar)
    debug && println("")
    debug && println("Performing the $(i)th MCMC iteration")
    if mod(i, 1) == 0
      step = rand(MvNormal(transition_cov))
    else
      step = fill(0., 6)
    end
    ilm_proposal = [ilm_trace.α[end],
                    ilm_trace.β[end],
                    ilm_trace.η[end],
                    ilm_trace.ρ[end],
                    ilm_trace.γ[end]] .+ step[1:5]
    detection_proposal = [detection_trace.ν[end]] .+ step[6]

    debug && println("Proposed parameter changes: $(round(step, 3))")
    debug && println("ILM proposal: $(round(ilm_proposal, 3))")
    debug && println("Detection proposal: $(round(detection_proposal, 3))")
    lp = logprior(ilm_priors, ilm_proposal, debug)
    lp += logprior(detection_priors, detection_proposal, debug)

    if lp > -Inf
      if mod(i, 1) == 0
        # Generate data augmentation proposal
        aug = propose_augment(ilm_proposal[4],
                              detection_proposal[1],
                              ilm_trace.network[end],
                              ilm_trace.aug[end],
                              obs,
                              rand(Poisson(1)),
                              debug)
      else
        aug = ilm_trace.aug[end]
      end
    end
    if lp > -Inf
      # SEIR loglikelihood
      ll, network_rates = SEIR_loglikelihood(ilm_proposal[1],
                                             ilm_proposal[2],
                                             ilm_proposal[3],
                                             ilm_proposal[4],
                                             ilm_proposal[5],
                                             aug,
                                             obs,
                                             debug,
                                             dist)
      lp += ll
    end
    if lp > -Inf
      if mod(i, 1) == 0
        # Generate network proposal
        network = propose_network(network_rates,
                                  ilm_trace.network[end],
                                  debug,
                                  rand(Poisson(1)),
                                  "multinomial")
      else
        network = ilm_trace.network[end]
      end
      lp += exposure_network_loglikelihood(network, network_rates, debug)
    end
    if lp > -Inf
      lp += detection_loglikelihood(detection_proposal[1],
                                   aug,
                                   network,
                                   obs,
                                   debug)
    end

    # Acceptance/rejection
    reject = MHreject(lp, ilm_trace.logposterior[end], debug)

    if reject
      aug = ilm_trace.aug[end]
      network_rates = ilm_trace.network_rates[end]
      network = ilm_trace.network[end]
      lp = ilm_trace.logposterior[end]
      ilm_proposal = [ilm_trace.α[end],
                      ilm_trace.β[end],
                      ilm_trace.η[end],
                      ilm_trace.ρ[end],
                      ilm_trace.γ[end]]
      detection_proposal = [detection_trace.ν[end]]
      rejects += 1
    end

    push!(ilm_trace.aug, aug)
    push!(ilm_trace.network_rates, network_rates)
    push!(ilm_trace.network, network)
    push!(ilm_trace.logposterior, lp)
    push!(ilm_trace.α, ilm_proposal[1])
    push!(ilm_trace.β, ilm_proposal[2])
    push!(ilm_trace.η, ilm_proposal[3])
    push!(ilm_trace.ρ, ilm_proposal[4])
    push!(ilm_trace.γ, ilm_proposal[5])
    push!(detection_trace.ν, detection_proposal[1])
  end
  println("MCMC acceptance rate: $(round(1.0-(rejects/n),4))")
  return ilm_trace, detection_trace
end


function MCMC(n::Int64,
              transition_cov::Array{Float64},
              ilm_priors::SEIR_priors,
              detection_priors::Lag_priors,
              obs::SEIR_observed,
              debug=false::Bool,
              progress=true::Bool,
              dist=Euclidean(),
              init_limit=500::Int64)

  ilm_trace, detection_trace = initialize(ilm_priors,
                                          detection_priors,
                                          obs,
                                          init_limit,
                                          debug,
                                          dist)
  return MCMC(n,
              transition_cov,
              ilm_trace,
              detection_trace,
              ilm_priors,
              detection_priors,
              obs,
              debug,
              progress,
              dist)
  end


function MCMC(n::Int64,
              ilm_priors::SEIR_priors,
              detection_priors::Lag_priors,
              obs::SEIR_observed,
              debug=false::Bool,
              progress=true::Bool,
              dist=Euclidean(),
              init_limit=500::Int64)

  transition_cov = diagm([var(ilm_priors.α),
                          var(ilm_priors.β),
                          var(ilm_priors.η),
                          var(ilm_priors.ρ),
                          var(ilm_priors.γ),
                          var(detection_priors.ν)])*(2.38^2)/6.
  return MCMC(n,
              transition_cov,
              ilm_priors,
              detection_priors,
              obs,
              debug,
              progress,
              dist,
              init_limit)
end


function MCMC(n::Int64,
              transition_cov::Array{Float64},
              ilm_trace::SEIR_trace,
              ilm_priors::SEIR_priors,
              obs::SEIR_observed,
              debug=false::Bool,
              progress=true::Bool,
              dist=Euclidean())

  @assert(size(transition_cov) == (5,5),
          "Transition kernel covariance matrix must be a positive definite 5x5 matrix")
  progressbar = Progress(n, 5, "Performing $n MCMC iterations...", 25)
  debug && println("MCMC transition kernel covariance matrix:")
  debug && println(round(transition_cov, 3))
  rejects = 0
  for i = 1:n
    progress && !debug && next!(progressbar)
    debug && println("")
    debug && println("Performing the $(i)th MCMC iteration")
    if mod(i, 1) == 0
      step = rand(MvNormal(transition_cov))
    else
      step = fill(0., 5)
    end
    ilm_proposal = [ilm_trace.α[end],
                    ilm_trace.β[end],
                    ilm_trace.η[end],
                    ilm_trace.ρ[end],
                    ilm_trace.γ[end]] .+ step[1:5]

    debug && println("Proposed parameter changes: $(round(step, 3))")
    debug && println("ILM proposal: $(round(ilm_proposal, 3))")
    lp = logprior(ilm_priors, ilm_proposal, debug)

    if lp > -Inf
      if mod(i, 2) == 1
        # Generate data augmentation proposal
        aug = propose_augment(ilm_proposal[4],
                              ilm_trace.network[end],
                              ilm_trace.aug[end],
                              obs,
                              rand(Poisson(2.)),
                              debug)
      else
        aug = ilm_trace.aug[end]
      end
    end
    if lp > -Inf
      # SEIR loglikelihood
      ll, network_rates = SEIR_loglikelihood(ilm_proposal[1],
                                             ilm_proposal[2],
                                             ilm_proposal[3],
                                             ilm_proposal[4],
                                             ilm_proposal[5],
                                             aug,
                                             obs,
                                             debug,
                                             dist)
      lp += ll
    end
    if lp > -Inf
      if mod(i, 2) == 0
        # Generate network proposal
        network = propose_network(network_rates,
                                  ilm_trace.network[end],
                                  debug,
                                  rand(Poisson(2.)),
                                  "uniform")
      else
        network = ilm_trace.network[end]
      end
      lp += exposure_network_loglikelihood(network, network_rates, debug)
    end

    # Acceptance/rejection
    reject = MHreject(lp, ilm_trace.logposterior[end], debug)

    if reject
      aug = ilm_trace.aug[end]
      network_rates = ilm_trace.network_rates[end]
      network = ilm_trace.network[end]
      lp = ilm_trace.logposterior[end]
      ilm_proposal = [ilm_trace.α[end],
                      ilm_trace.β[end],
                      ilm_trace.η[end],
                      ilm_trace.ρ[end],
                      ilm_trace.γ[end]]
      rejects += 1
    end

    push!(ilm_trace.aug, aug)
    push!(ilm_trace.network_rates, network_rates)
    push!(ilm_trace.network, network)
    push!(ilm_trace.logposterior, lp)
    push!(ilm_trace.α, ilm_proposal[1])
    push!(ilm_trace.β, ilm_proposal[2])
    push!(ilm_trace.η, ilm_proposal[3])
    push!(ilm_trace.ρ, ilm_proposal[4])
    push!(ilm_trace.γ, ilm_proposal[5])
  end
  println("MCMC acceptance rate: $(round(1.0-(rejects/n),4))")
  return ilm_trace
end


function MCMC(n::Int64,
              transition_cov::Array{Float64},
              ilm_priors::SEIR_priors,
              obs::SEIR_observed,
              debug=false::Bool,
              progress=true::Bool,
              dist=Euclidean(),
              init_limit=500::Int64)

  ilm_trace = initialize(ilm_priors,
                         obs,
                         init_limit,
                         debug,
                         dist)
  return MCMC(n,
              transition_cov,
              ilm_trace,
              ilm_priors,
              obs,
              debug,
              progress,
              dist)
  end


function MCMC(n::Int64,
              ilm_priors::SEIR_priors,
              obs::SEIR_observed,
              debug=false::Bool,
              progress=true::Bool,
              dist=Euclidean(),
              init_limit=500::Int64)

  transition_cov = diagm([var(ilm_priors.α),
                          var(ilm_priors.β),
                          var(ilm_priors.η),
                          var(ilm_priors.ρ),
                          var(ilm_priors.γ)])*(2.38^2)/5.
  return MCMC(n,
              transition_cov,
              ilm_priors,
              obs,
              debug,
              progress,
              dist,
              init_limit)
end
