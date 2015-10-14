"""
infer.jl
"""

function surveil(population::Population, ν::Float64)
  """
  Gather surveillance data on specific individuals in a population, with an exponentially distributed detection lag with rate ν
  """
  exposed_actual = fill(NaN, length(population.events)-1)
  infectious_actual = fill(NaN, length(population.events)-1)
  infectious_observed = fill(NaN, length(population.events)-1)
  removed_actual = fill(NaN, length(population.events)-1)
  removed_observed = fill(NaN, length(population.events)-1)
  covariates_actual = fill(fill(NaN, length(population.history[2][1][1])),  length(population.events)-1)
  covariates_observed = fill(fill(NaN, length(population.history[2][1][1])),  length(population.events)-1)
  seq_actual = convert(Vector{Any}, fill(NaN, length(population.events)-1))
  seq_observed = convert(Vector{Any}, fill(NaN, length(population.events)-1))

  for i = 2:length(population.events)
    # Initial conditions (assumed to be constant and observed without error)
    covariates_actual[i-1] = population.history[i][1][1]
    covariates_observed[i-1] = population.history[i][1][1]

    # Exposure time (unobservable)
    if length(population.events[i][1]) > 0
      exposed_actual[i-1] = population.events[i][1][1]
    end

    # Infectious time (observed with latency)
    if length(population.events[i][3]) > 0
      infectious_actual[i-1] = population.events[i][3][1]
      seq_actual[i-1] = convert(Vector{Int64}, population.history[i][2][find(infectious_actual[i-1] .>= population.events[i][6])[end]])
      if ν < Inf
        infectious_observed[i-1] = infectious_actual[i-1] + rand(Exponential(1/ν))
      elseif ν == Inf
        infectious_observed[i-1] = infectious_actual[i-1]
      end

      if length(population.events[i][4]) > 0 && infectious_observed[i-1] >= population.events[i][4][1]
        infectious_observed[i-1] = NaN
      else
        seq_observed[i-1] = convert(Vector{Int64}, population.history[i][2][find(infectious_observed[i-1] .>= population.events[i][6])[end]])
      end
    end

    # Removal time (observed with latency)
    if length(population.events[i][4]) > 0
      removed_actual[i-1] = population.events[i][4][1]
      if !isnan(infectious_observed[i-1])
        if ν < Inf
          removed_observed[i-1] = removed_actual[i-1] + rand(Exponential(1/ν))
        elseif ν == Inf
          removed_observed[i-1] = removed_actual[i-1]
        end
      end
    end
  end
  return SEIR_actual(exposed_actual, infectious_actual, removed_actual, covariates_actual, seq_actual), SEIR_observed(infectious_observed, removed_observed, covariates_observed, seq_observed)
end


surveil(population, Inf) = surveil(population::Population)


function propose_augment(changed_individuals::Vector{Int64}, ρ::Float64, ν::Float64, network::Array{Bool, 2}, previous_aug::SEIR_augmented, obs::SEIR_observed, debug=false::Bool)
  """
  Proposes augmented data for a specified vector of `changed_individuals`
  """
  exposed_augmented = copy(previous_aug.exposed)
  infectious_augmented = copy(previous_aug.infectious)
  removed_augmented = copy(previous_aug.removed)

  for i = changed_individuals
    if ν < Inf
      pathway_out = pathwayfrom(i, network)
      pathway_in = pathwayto(i, network)
      if length(pathway_in) > 2
        infectious_augmented[i] = obs.infectious[i] - rand(Truncated(Exponential(1/ν), obs.infectious[i] - minimum(obs.infectious[pathway_out]), obs.infectious[i] - infectious_augmented[pathway_in[2]]))
        if isnan(obs.removed[pathway_in[2]])
          exposed_augmented[i] = infectious_augmented[i] - rand(Truncated(Exponential(1/ρ), 0., infectious_augmented[i]-infectious_augmented[pathway_in[2]]))
        else
          exposed_augmented[i] = infectious_augmented[i] - rand(Truncated(Exponential(1/ρ), infectious_augmented[i]-removed_augmented[pathway_in[2]], infectious_augmented[i]-infectious_augmented[pathway_in[2]]))
        end
      else
        infectious_augmented[i] = obs.infectious[i] - rand(Truncated(Exponential(1/ν), obs.infectious[i] - minimum(obs.infectious[pathway_out]), Inf))
        exposed_augmented[i] = infectious_augmented[i] - rand(Exponential(1/ρ))
      end
      if !isnan(obs.removed[i])
        removed_augmented[i] = obs.removed[i] - rand(Truncated(Exponential(1/ν), 0., obs.removed[i] - obs.infectious[i]))
      end
    elseif ν == Inf
      infectious_augmented[i] = obs.infectious[i]
      if isnan(obs.removed[pathway_in[2]])
        exposed_augmented[i] = infectious_augmented[i] - rand(Truncated(Exponential(1/ρ), 0., infectious_augmented[i]-infectious_augmented[pathway_in[2]]))
      else
        exposed_augmented[i] = infectious_augmented[i] - rand(Truncated(Exponential(1/ρ), infectious_augmented[i]-removed_augmented[pathway_in[2]], infectious_augmented[i]-infectious_augmented[pathway_in[2]]))
      end
      if !isnan(obs.removed[i])
        removed_augmented[i] = obs.removed[i]
      end
    end
  end
  if debug
    println("DATA AUGMENTATION")
    println("$(sum(!isnan(exposed_augmented))), $(sum(!isnan(infectious_augmented))), and $(sum(!isnan(removed_augmented))) augmented exposure, infection, and removal times respectively")
    for i = 1:length(obs.infectious)
      @assert(!(isnan(exposed_augmented[i]) && !isnan(obs.infectious[i])), "Data augmentation error: could not generate exposure event $i")
      @assert(!(isnan(infectious_augmented[i]) && !isnan(obs.infectious[i])), "Data augmentation error: could not generate infectious event $i")
      @assert(!(isnan(removed_augmented[i]) && !isnan(obs.removed[i])), "Data augmentation error: could not generate exposure event $i")
    end
  end
  return SEIR_augmented(exposed_augmented, infectious_augmented, removed_augmented)
end


function propose_augment(ρ::Float64, ν::Float64, network::Array{Bool, 2}, previous_aug::SEIR_augmented, obs::SEIR_observed, changes=1::Int64, debug=false::Bool)
  """
  Proposes augmented data by making random selection of individual `changes` to augmented data
  """
  if changes == 0
    changed_individuals = pathwayfrom(0, network)
  else
    changed_individuals = sample(pathwayfrom(0, network), changes, replace=false)
  end
  return propose_augment(changed_individuals, ρ, ν, network, previous_aug, obs, debug)
end

propose_augment(ρ, Inf, network, previous_aug, obs, changes, debug) = propose_augment(ρ::Float64, network::Array{Bool, 2}, previous_aug::SEIR_augmented, obs::SEIR_observed, changes=1::Int64, debug=false::Bool)

function propose_augment(ρ::Float64, ν::Float64, obs::SEIR_observed, debug=false::Bool)
  """
  Augments surveilance data, organizes observations
  """
  exposed_augmented = fill(NaN, length(obs.infectious))
  infectious_augmented = fill(NaN, length(obs.infectious))
  removed_augmented = fill(NaN, length(obs.removed))
  for i = 1:length(obs.infectious)
    if !isnan(obs.infectious[i])
      if ν < Inf
        infectious_augmented[i] = obs.infectious[i] - rand(Exponential(1/ν))
      elseif ν == Inf
        infectious_augmented[i] = obs.infectious[i]
      end
      exposed_augmented[i] = infectious_augmented[i] - rand(Exponential(1/ρ))
      if !isnan(obs.removed[i])
        if ν < Inf
          removed_augmented[i] = obs.removed[i] - rand(Truncated(Exponential(1/ν), -Inf, obs.removed[i] - obs.infectious[i]))
        elseif ν == Inf
          removed_augmented[i] = obs.removed[i]
        end
      end
    end
  end
  if debug
    println("DATA AUGMENTATION")
    println("Augmented exposure times ($(sum(!isnan(exposed_augmented))) total): $(round(exposed_augmented,3))")
    println("Augmented infection times ($(sum(!isnan(removed_augmented))) total): $(round(infectious_augmented,3))")
    println("Augmented removal times ($(sum(!isnan(removed_augmented))) total): $(round(removed_augmented,3))")
    println("")
  end
  return SEIR_augmented(exposed_augmented, infectious_augmented, removed_augmented)
end


propose_augment(ρ, Inf, obs, debug) = propose_augment(ρ::Float64, obs::SEIR_observed, debug=false::Bool)


function logprior(priors::Priors, params::Vector{Float64})
  """
  Calculate the logprior from prior distributions
  """
  lprior = 0.
  for i = 1:length(params)
    lprior += logpdf(priors.(names(priors)[i]), params[i])
  end
  return lprior
end


function randprior(priors::Priors)
  """
  Randomly generate a parameter vector from specified priors
  """
  params = [rand(priors.(names(priors)[1]))]
  for i = 2:length(names(priors))
    push!(params, rand(priors.(names(priors)[i])))
  end
  return params
end


function propose_network(network_rates::Array{Float64, 2}, previous_network::Array{Bool, 2}, debug=false::Bool, changes=1::Int64, method="multinomial"::String)
  """
  Propose a network
  """
  @assert(changes <= sum(rate_totals .> 0), "Attempting to make more network changes than there are exposure events")
  if changes == 0
    changed_individuals = find(rate_totals .> 0)
  else
    changed_individuals = sample(find(rate_totals .> 0), changes, replace=false)
  end
  return propose_network(changed_individuals, network_rates, previous_network, debug, method)
end


function propose_network(changed_individuals::Vector{Int64}, network_rates::Array{Float64, 2}, previous_network::Array{Bool, 2}, debug=false::Bool, method="multinomial"::String)
  """
  Propose a network
  """
  @assert(any(method .== ["uniform", "multinomial"]), "Network proposal method must be 'uniform' or 'multinomial'.")
  network = copy(previous_network)
  rate_totals = sum(network_rates,1)
  network[:, changed_individuals] = false
  if method == "uniform"
    for i = changed_individuals
      network[sample(find(network_rates[:,i] .> 0.)), i] = true
    end
  elseif method == "multinomial"
    @assert(size(network_rates) == size(previous_network), "A mismatch in the previous network and network rates dimensions was detected in the network proposal function")
    for i = changed_individuals
      network[findfirst(rand(Multinomial(1, network_rates[:,i]/rate_totals[i]))), i] = true
    end
  end
  if debug
    println("Network proposal ($(sum(network)) infections total, up to $(length(changed_individuals)) changes from previous network):")
    println("$(0 + network)")
  end
  return network
end


function propose_network(network_rates::Array{Float64, 2}, debug=false::Bool)
  """
  Initial network proposal
  """
  network = fill(false, size(network_rates))
  rate_totals = sum(network_rates,1)
  exposures = find(rate_totals .> 0)
  for i = exposures
    network[findfirst(rand(Multinomial(1, network_rates[:,i]/rate_totals[i]))), i] = true
  end
  if debug
    println("Network proposal ($(sum(network)) infections total):")
    println("$(0 + network)")
  end
  return network
end


function seq_distances(obs::SEIR_observed, aug::SEIR_augmented, network::Array{Bool, 2}, debug=false::Bool)
  """
  For a given transmission network, find the time between the pathogen sequences between every individuals i and j
  """
  pathways = pathwaysto(network)

  seq_dist = fill(0., (size(network, 2), size(network, 2)))

  for i = 1:length(pathways)
    if debug
      println("")
      println("SEQUENCE DISTANCES")
      println("Infection of individual $(pathways[i][1]) observed at $(obs.infectious[pathways[i][1]])")
      println("Infection pathway of individual $(pathways[i][1]) is $(pathways[i])")
    end
    for j = 1:(i-1)
      k = 1
      while length(pathways[i]) > k && length(pathways[j]) > k && pathways[i][end - k] == pathways[j][end - k]
        k += 1
      end
      if debug
        println("Infection of individual $(pathways[j][1]) observed at $(obs.infectious[pathways[j][1]])")
        println("Infection pathway of individual $(pathways[j][1]) is $(pathways[j])")
        if k == length(pathways[i]) || k == length(pathways[j])
          println("Linear infection pathway between individual $(pathways[i][1]) and individual $(pathways[j][1])")
        else
          println("Most recent common infection source of individuals $(pathways[i][1]) and $(pathways[j][1]) is $(pathways[i][end - k + 1])")
          println("The infection pathway of $(pathways[i][1]) and $(pathways[j][1]) diverged with $(pathways[i][end - k]) and $(pathways[j][end - k])")
        end
      end

      if k == length(pathways[i])
        seq_dist[pathways[i][1],pathways[j][1]] += obs.infectious[pathways[j][1]] - aug.exposed[pathways[j][end - k]]
        seq_dist[pathways[i][1],pathways[j][1]] += abs(aug.exposed[pathways[j][end - k]] - obs.infectious[pathways[i][1]])
      elseif k == length(pathways[j])
        seq_dist[pathways[i][1],pathways[j][1]] += obs.infectious[pathways[i][1]] - aug.exposed[pathways[i][end - k]]
        seq_dist[pathways[i][1],pathways[j][1]] += abs(aug.exposed[pathways[i][end - k]] - obs.infectious[pathways[j][1]])
      else
        seq_dist[pathways[i][1],pathways[j][1]] += obs.infectious[pathways[i][1]] - aug.exposed[pathways[i][end - k]]
        seq_dist[pathways[i][1],pathways[j][1]] += obs.infectious[pathways[j][1]] - aug.exposed[pathways[j][end - k]]
        seq_dist[pathways[i][1],pathways[j][1]] += abs(aug.exposed[pathways[j][end - k]] - aug.exposed[pathways[i][end - k]])
      end
    end
    if debug
      println("")
    end
  end

  return seq_dist += transpose(seq_dist)
end


function network_loglikelihood(obs::SEIR_observed, aug::SEIR_augmented, network::Array{Bool, 2}, p_matrix::Function, debug=false::Bool)
  """
  Loglikelihood for an entire transmission network
  """
  ll = 0.
  infected = find(sum(network, 1))
  seq_dist = seq_distances(obs, aug, network, debug)

  for i = 1:length(infected)
    for j = 1:(i-1)
      ll += sum(log(p_matrix(seq_dist[infected[i],infected[j]]))[sub2ind((4,4), obs.seq[infected[i]], obs.seq[infected[j]])])
    end
  end
  return ll
end


function SEIR_loglikelihood(α::Float64, β::Float64, η::Float64, ρ::Float64, γ::Float64, aug::SEIR_augmented, obs::SEIR_observed, debug=false::Bool, dist=Euclidean())
  """
  Calculate the loglikelihood and return an exposure network array under specified parameters values and observations

  α, β: powerlaw exposure kernel parameters
  η: external pressure rate
  ρ: infectivity rate (1/mean latent period)
  γ: removal rate (1/mean infectious period)
  """
  # Initiate an exposure network
  network_rates = fill(0., (1 + length(obs.covariates), length(obs.covariates)))

  # Start with loglikelihood of 0
  ll = 0.

  # Create event timing array
  event_times = [aug.exposed aug.infectious aug.removed]

  # Find event order
  event_order = sortperm(event_times[:])

  # Create empty rate array
  rate_array = fill(0., (1 + length(obs.covariates) + 2, length(obs.covariates)))

  # First row is external pressure rate
  rate_array[1,:] = η

  # Loglikelihood of events that have been observed
  for i = 1:length(event_order)

    # Stop loglikelihood calculation after last event considered
    isnan(event_times[event_order[i]]) && break

    # Stop loglikelihood calculation anytime the loglikelihood goes to -Inf
    ll == -Inf && break

    # Convert linear index to an event tuple (individual, event type)
    id = ind2sub(size(event_times), event_order[i])

    # Don't consider likelilihood contribution of first event
    if i > 1
      # Find the "master rate" through the sum of rate_array
      master_rate = sum(rate_array)

      # loglikelihood of event time with master rate
      ll += logpdf(Exponential(1/master_rate), event_times[event_order[i]] - event_times[event_order[i-1]])

      # loglikelihood of the one event that did occur
      ll += log(sum(rate_array[:, id[1]]) / master_rate)
    end

    # Exposure event
    if id[2] == 1

      # Record exposure rates at time of exposure
      network_rates[:, id[1]] = copy(rate_array[1:(length(obs.covariates)+1), id[1]])

      # Update exposure rates
      rate_array[1:(1 + size(rate_array, 2)), id[1]] = 0.

      # Update infectivity rate
      rate_array[1 + size(rate_array, 2) + 1, id[1]] = ρ

    # Infectiousness event
    elseif id[2] == 2
      # Update infectivity rate
      rate_array[1 + size(rate_array, 2) + 1, id[1]] = 0.

      # Update removal rate
      rate_array[1 + size(rate_array, 2) + 2, id[1]] = γ

      # Update exposure rates for rest of susceptible population
      for j = 1:size(rate_array, 2)
        if j != id[1] && rate_array[1, j] > 0.
          rate_array[id[1] + 1, j] = α*evaluate(dist, obs.covariates[id[1]], obs.covariates[j])^-β
        end
      end

    # Removal event
    elseif id[2] == 3
      # Update removal rate
      rate_array[1 + size(rate_array, 2) + 2, id[1]] = 0.

      # Update exposure rates
      rate_array[id[1] + 1, :] = 0.
    end

    # Provide loop position when loglikelihood goes to -Inf when debugging
    if debug && ll == -Inf
      println("SEIR LOGLIKELIHOOD")
      if id[2] == 1
        println("Event $i (exposure of individual $(id[1])) caused loglikelihood to go to -Inf")
      elseif id[2] == 2
        println("Event $i (infection of individual $(id[1])) caused loglikelihood to go to -Inf")
      elseif id[2] == 3
        println("Event $i (removal of individual $(id[1])) caused loglikelihood to go to -Inf")
      end
      println("")
    end
  end
  return ll, network_rates
end


function initialize(ilm_priors::SEIR_priors, mutation_priors::JC69_priors, detection_priors::Lag_priors, obs::SEIR_observed, limit=500::Int, debug=false::Bool, dist=Euclidean())
  """
  Initiate an Trace object by sampling from specified prior distributions
  """
  ilm_params = randprior(ilm_priors)
  mutation_params = randprior(mutation_priors)
  detection_params = randprior(detection_priors)
  aug = augment(ilm_params[4], detection_params[1], obs, debug)
  lp1, network_rates = SEIR_loglikelihood(ilm_params[1], ilm_params[2], ilm_params[3], ilm_params[4], ilm_params[5], aug, obs, debug, dist)
  count = 1

  # Retry initialization until non-negative infinity loglikelihood
  while lp1 == -Inf && count < limit
    count += 1
    ilm_params = randprior(ilm_priors)
    mutation_params = randprior(mutation_priors)
    detection_params = randprior(detection_priors)
    aug = augment(ilm_params[4], detection_params[1], obs, debug)
    lp1, network_rates = SEIR_loglikelihood(ilm_params[1], ilm_params[2], ilm_params[3], ilm_params[4], ilm_params[5], aug, obs, debug, dist)
    lp1 += logprior(ilm_priors, ilm_params) + logprior(detection_priors, detection_params)
  end

  if count < limit
    network = propose_network(network_rates, debug)
    lp2 = network_loglikelihood(obs, aug, network, jc69p([mutation_params[1]]), debug)
    lp2 += logprior(mutation_priors, mutation_params)

    while lp1 + lp2 == -Inf && count < limit
      count += 1
      lp1, network_rates = SEIR_loglikelihood(ilm_params[1], ilm_params[2], ilm_params[3], ilm_params[4], ilm_params[5], aug, obs, debug, dist)
      lp1 += logprior(ilm_priors, ilm_params) + logprior(detection_priors, detection_params)
      network = propose_network(network_rates, debug)
      lp2 = network_loglikelihood(obs, aug, network, jc69p([mutation_params[1]]), debug)
      lp2 += logprior(mutation_priors, mutation_params)
    end

    if count < limit && lp1 + lp2 > -Inf
      println("INITIALIZATON")
      println("Successful on attempt $count (lp1 = $lp1, lp2 = $lp2)")
      return SEIR_trace([ilm_params[1]], [ilm_params[2]], [ilm_params[3]], [ilm_params[4]], [ilm_params[5]], [aug], Array[network_rates], Array[network], [lp1], [lp2]), Lag_trace([detection_params[1]]), JC69_trace([mutation_params[1]])
    else
      println("INITIALIZATON")
      println("Failed to initialize after $count attempts (lp1 = $lp1, lp2 = $lp2)")
    end
  end
end


function initialize(ilm_priors::SEIR_priors, detection_priors::Lag_priors, obs::SEIR_observed, limit=500::Int, debug=false::Bool, dist=Euclidean())
  """
  Initiate an Trace object by sampling from specified prior distributions
  """
  ilm_params = randprior(ilm_priors)
  detection_params = randprior(detection_priors)
  aug = augment(ilm_params[4], detection_params[1], obs, debug)
  lp1, network_rates = SEIR_loglikelihood(ilm_params[1], ilm_params[2], ilm_params[3], ilm_params[4], ilm_params[5], aug, obs, debug, dist)
  lp1 += logprior(ilm_priors, ilm_params) + logprior(detection_priors, detection_params)
  count = 1

  # Retry initialization until non-negative infinity loglikelihood
  while lp1 == -Inf && count < limit
    count += 1
    ilm_params = randprior(ilm_priors)
    detection_params = randprior(detection_priors)
    aug = augment(ilm_params[4], detection_params[1], obs, debug)
    lp1, network_rates = SEIR_loglikelihood(ilm_params[1], ilm_params[2], ilm_params[3], ilm_params[4], ilm_params[5], aug, obs, debug, dist)
    lp1 += logprior(ilm_priors, ilm_params) + logprior(detection_priors, detection_params)
  end

  if count < limit
    network = propose_network(network_rates, debug)
    println("Successful initalization on attempt $count (marginal log posterior = $lp1)")
    return SEIR_trace([ilm_params[1]], [ilm_params[2]], [ilm_params[3]], [ilm_params[4]], [ilm_params[5]], [aug], Array[network_rates], Array[network], [lp1], [0]), Lag_trace([detection_params[1]])
  else
    println("Failed to initialize after $count attempts (marginal log posterior = $lp1)")
  end
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
  """
  Performs `n` data-augmented metropolis hastings within Gibbs MCMC iterations. Initiates a single chain by sampling from prior distribution
  """

  @assert(size(transition_cov) == (7,7), "Transition kernel's covariance matrix must be a positive definite 7x7 matrix")

  for i = 1:n
    # Create and incremenet progress bar
    if progress
      if i == 1
        progressbar = Progress(n, 5, "Performing $n MCMC iterations...", 30)
      else
        next!(progressbar)
      end
    end

    # Step 2a: Metropolis-Hastings proposal
    # Only generate valid proposals
    step = rand(MvNormal(transition_cov))
    ilm_proposal = [ilm_trace.α[end], ilm_trace.β[end], ilm_trace.η[end], ilm_trace.ρ[end], ilm_trace.γ[end]] .+ step[1:5]
    detection_proposal = [detection_trace.ν[end]] .+ step[6]
    mutation_proposal = [mutation_trace.λ[end]] .+ step[7]

    lp1b_proposal = logprior(ilm_priors, ilm_proposal)
    lp1b_proposal += logprior(detection_priors, detection_proposal)
    lp2b_proposal = logprior(mutation_priors, mutation_proposal)

    while lp1b_proposal + lp2b_proposal == -Inf
      step = rand(MvNormal(transition_cov))
      ilm_proposal = [ilm_trace.α[end], ilm_trace.β[end], ilm_trace.η[end], ilm_trace.ρ[end], ilm_trace.γ[end]] .+ step[1:5]
      detection_proposal = [detection_trace.ν[end]] .+ step[6]
      mutation_proposal = [mutation_trace.λ[end]] .+ step[7]
      lp1b_proposal = logprior(ilm_priors, ilm_proposal)
      lp1b_proposal += logprior(detection_priors, detection_proposal)
      lp2b_proposal = logprior(mutation_priors, mutation_proposal)
    end

    # Propose new augmented data
    aug_proposal = augment(ilm_trace.ρ[end], detection_proposal[1], ilm_trace.network[end], obs, debug)
    lp1a_proposal, network_rates_proposal = SEIR_loglikelihood(ilm_trace.α[end], ilm_trace.β[end], ilm_trace.η[end], ilm_trace.ρ[end], ilm_trace.γ[end], aug_proposal, obs, debug, dist)
    lp1a_proposal += logprior(ilm_priors, [ilm_trace.α[end], ilm_trace.β[end], ilm_trace.η[end], ilm_trace.ρ[end], ilm_trace.γ[end]])
    lp1a_proposal += logprior(detection_priors, [detection_proposal])
    lp2a_proposal = network_loglikelihood(obs, aug_proposal, ilm_trace.network[end], jc69p([mutation_trace.λ[end]]), debug)
    lp2a_proposal += logprior(mutation_priors, [mutation_trace.λ[end]])

    reject = true
    if lp1a_proposal + lp2a_proposal >= ilm_trace.logposterior_1[end] + ilm_trace.logposterior_2[end]
      reject = false
    elseif exp(lp1a_proposal + lp2a_proposal - ilm_trace.logposterior_1[end] + ilm_trace.logposterior_2[end]) >= rand()
      reject = false
    end

    if reject
      aug_proposal = ilm_trace.aug[end]
      network_rates_proposal = ilm_trace.network_rates[end]
      lp1a_proposal = ilm_trace.logposterior_1[end]
      lp2a_proposal = ilm_trace.logposterior_2[end]
      detection_proposal = detection_trace.ν[end]
    end

    push!(ilm_trace.aug, aug_proposal)
    push!(ilm_trace.network_rates, network_rates_proposal)
    push!(ilm_trace.logposterior_1, lp1a_proposal)
    push!(ilm_trace.logposterior_2, lp2a_proposal)
    push!(detection_trace.ν, detection_proposal[1])

    # Propose new disease model parameters
    ll_proposal, network_rates_proposal = SEIR_loglikelihood(ilm_proposal[1], ilm_proposal[2], ilm_proposal[3], ilm_proposal[4], ilm_proposal[5], ilm_trace.aug[end], obs, debug, dist)
    lp1b_proposal += ll_proposal

    reject = true
    if lp1b_proposal >= ilm_trace.logposterior_1[end]
      reject = false
    elseif exp(lp1b_proposal - ilm_trace.logposterior_1[end]) >= rand()
      reject = false
    end

    if reject
      ilm_proposal[1] = ilm_trace.α[end]
      ilm_proposal[2] = ilm_trace.β[end]
      ilm_proposal[3] = ilm_trace.η[end]
      ilm_proposal[4] = ilm_trace.ρ[end]
      ilm_proposal[5] = ilm_trace.γ[end]
      network_rates_proposal = ilm_trace.network_rates[end]
      lp1b_proposal = ilm_trace.logposterior_1[end]
    end

    push!(ilm_trace.α, ilm_proposal[1])
    push!(ilm_trace.β, ilm_proposal[2])
    push!(ilm_trace.η, ilm_proposal[3])
    push!(ilm_trace.ρ, ilm_proposal[4])
    push!(ilm_trace.γ, ilm_proposal[5])
    ilm_trace.network_rates[end] = network_rates_proposal
    ilm_trace.logposterior_1[end] = lp1b_proposal

    # Propose new network
    network_proposal = propose_network(ilm_trace.network_rates[end], ilm_trace.network[end], debug)
    lp2b_proposal += network_loglikelihood(obs, aug_proposal, network_proposal, jc69p([mutation_proposal]), debug)

    reject = true
    if lp2b_proposal >= ilm_trace.logposterior_2[end]
      reject = false
    elseif exp(lp2b_proposal - ilm_trace.logposterior_2[end]) >= rand()
      reject = false
    end

    if reject
      network_proposal = ilm_trace.network[end]
      lp2b_proposal = ilm_trace.logposterior_2[end]
      mutation_proposal = [mutation_trace.λ[end]]
    end

    push!(ilm_trace.network, network_proposal)
    ilm_trace.logposterior_2[end] = lp2b_proposal
    push!(mutation_trace.λ, mutation_proposal[1])
  end
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
  """
  Performs `n` data-augmented metropolis hastings within Gibbs MCMC iterations. Initiates a single chain by sampling from prior distribution
  """

  ilm_trace, detection_trace, mutation_trace = initialize(ilm_priors, mutation_priors, detection_priors, obs, init_limit, debug, dist)

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
  """
  Performs `n` data-augmented metropolis hastings within Gibbs MCMC iterations. Initiates a single chain by sampling from prior distribution
  """

  transition_cov = diagm([var(ilm_priors.α),
                          var(ilm_priors.β),
                          var(ilm_priors.η),
                          var(ilm_priors.ρ),
                          var(ilm_priors.γ),
                          var(detection_priors.ν),
                          var(mutation_priors.λ)])

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
  """
  Performs `n` data-augmented metropolis hastings within Gibbs MCMC iterations. Initiates a single chain by sampling from prior distribution
  """

  @assert(size(transition_cov) == (6,6), "Transition kernel's covariance matrix must be a positive definite 6x6 matrix")

  for i = 1:n

    # Create and incremenet progress bar
    if progress
      if i == 1
        progressbar = Progress(n, 5, "Performing $n MCMC iterations...", 30)
      else
        next!(progressbar)
      end
    end

    # Only generate valid Metropolis-Hastings proposals
    step = rand(MvNormal(transition_cov))
    ilm_proposal = [ilm_trace.α[end], ilm_trace.β[end], ilm_trace.η[end], ilm_trace.ρ[end], ilm_trace.γ[end]] .+ step[1:5]
    detection_proposal = [detection_trace.ν[end]] .+ step[6]

    lp1b_proposal = logprior(ilm_priors, ilm_proposal)
    lp1b_proposal += logprior(detection_priors, detection_proposal)

    while lp1b_proposal == -Inf
      step = rand(MvNormal(transition_cov))
      ilm_proposal = [ilm_trace.α[end], ilm_trace.β[end], ilm_trace.η[end], ilm_trace.ρ[end], ilm_trace.γ[end]] .+ step[1:5]
      detection_proposal = [detection_trace.ν[end]] .+ step[6]
      lp1b_proposal = logprior(ilm_priors, ilm_proposal)
      lp1b_proposal += logprior(detection_priors, detection_proposal)
    end

    # Propose new augmented data
    aug_proposal = augment(ilm_trace.ρ[end], detection_proposal[1], ilm_trace.network[end], obs, debug)
    lp1a_proposal, network_rates_proposal = SEIR_loglikelihood(ilm_trace.α[end], ilm_trace.β[end], ilm_trace.η[end], ilm_trace.ρ[end], ilm_trace.γ[end], aug_proposal, obs, debug, dist)
    lp1a_proposal += logprior(ilm_priors, [ilm_trace.α[end], ilm_trace.β[end], ilm_trace.η[end], ilm_trace.ρ[end], ilm_trace.γ[end]])
    lp1a_proposal += logprior(detection_priors, [detection_proposal])

    reject = true
    if lp1a_proposal >= ilm_trace.logposterior_1[end]
      reject = false
    elseif exp(lp1a_proposal - ilm_trace.logposterior_1[end]) >= rand()
      reject = false
    end

    if reject
      aug_proposal = ilm_trace.aug[end]
      network_rates_proposal = ilm_trace.network_rates[end]
      lp1a_proposal = ilm_trace.logposterior_1[end]
      detection_proposal = detection_trace.ν[end]
    end

    push!(ilm_trace.aug, aug_proposal)
    push!(ilm_trace.network_rates, network_rates_proposal)
    push!(ilm_trace.logposterior_1, lp1a_proposal)
    push!(detection_trace.ν, detection_proposal[1])

    # Propose new disease model parameters
    ll_proposal, network_rates_proposal = SEIR_loglikelihood(ilm_proposal[1], ilm_proposal[2], ilm_proposal[3], ilm_proposal[4], ilm_proposal[5], ilm_trace.aug[end], obs, debug, dist)
    lp1b_proposal += ll_proposal

    reject = true
    if lp1b_proposal >= ilm_trace.logposterior_1[end]
      reject = false
    elseif exp(lp1b_proposal - ilm_trace.logposterior_1[end]) >= rand()
      reject = false
    end

    if reject
      ilm_proposal[1] = ilm_trace.α[end]
      ilm_proposal[2] = ilm_trace.β[end]
      ilm_proposal[3] = ilm_trace.η[end]
      ilm_proposal[4] = ilm_trace.ρ[end]
      ilm_proposal[5] = ilm_trace.γ[end]
      network_rates_proposal = ilm_trace.network_rates[end]
      lp1b_proposal = ilm_trace.logposterior_1[end]
    end

    push!(ilm_trace.α, ilm_proposal[1])
    push!(ilm_trace.β, ilm_proposal[2])
    push!(ilm_trace.η, ilm_proposal[3])
    push!(ilm_trace.ρ, ilm_proposal[4])
    push!(ilm_trace.γ, ilm_proposal[5])
    ilm_trace.network_rates[end] = network_rates_proposal
    ilm_trace.logposterior_1[end] = lp1b_proposal

    # Propose new network
    network_proposal = propose_network(ilm_trace.network_rates[end], ilm_trace.network[end], debug)
    lp2_proposal = 0.

    reject = true
    if lp2_proposal >= ilm_trace.logposterior_2[end]
      reject = false
    elseif exp(lp2_proposal - ilm_trace.logposterior_2[end]) >= rand()
      reject = false
    end

    if reject
      network_proposal = ilm_trace.network[end]
      lp2_proposal = ilm_trace.logposterior_2[end]
    end

    push!(ilm_trace.network, network_proposal)
    push!(ilm_trace.logposterior_2, lp2_proposal)
  end

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
  """
  Performs `n` data-augmented metropolis hastings within Gibbs MCMC iterations. Initiates a single chain by sampling from prior distribution
  """

  ilm_trace, detection_trace = initialize(ilm_priors, detection_priors, obs, init_limit, debug, dist)

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
  """
  Performs `n` data-augmented metropolis hastings within Gibbs MCMC iterations. Initiates a single chain by sampling from prior distribution
  """

  transition_cov = diagm([var(ilm_priors.α),
                          var(ilm_priors.β),
                          var(ilm_priors.η),
                          var(ilm_priors.ρ),
                          var(ilm_priors.γ),
                          var(detection_priors.ν)])

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
