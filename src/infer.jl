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
      seq_actual[i-1] = population.history[i][2][find(infectious_actual[i-1] .>= population.events[i][6])[end]]
      if ν < Inf
        infectious_observed[i-1] = infectious_actual[i-1] + rand(Exponential(1/ν))
      elseif ν == Inf
        infectious_observed[i-1] = infectious_actual[i-1]
      end

      if length(population.events[i][4]) > 0 && infectious_observed[i-1] >= population.events[i][4][1]
        infectious_observed[i-1] = NaN
      else
        seq_observed[i-1] = population.history[i][2][find(infectious_observed[i-1] .>= population.events[i][6])[end]]
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

surveil(population::Population) = surveil(population, Inf)

function augment(ρ::Float64, ν::Float64, obs::SEIR_observed)
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
  return SEIR_augmented(exposed_augmented, infectious_augmented, removed_augmented)
end

augment(ρ::Float64, obs::SEIR_observed) = augment(ρ, Inf, obs)

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

function seq_distances(obs::SEIR_observed, aug::SEIR_augmented, network::Array, debug=false::Bool)
  """
  For a given transmission network, find the time between the pathogen sequences between every individuals i and j
  """
  infected = find(isseq(obs.seq))
  pathway = [infected[1]]
  while pathway[end] != 0
    push!(pathway, findfirst(network[:,pathway[end]])-1)
  end
  pathways = Vector[pathway]

  for i = 2:length(infected)
    pathway = [infected[i]]
    while pathway[end] != 0
      push!(pathway, findfirst(network[:,pathway[end]])-1)
    end
    push!(pathways, pathway)
  end

  seq_dist = fill(0., (size(network, 2), size(network, 2)))

  for i = 1:length(infected)
    for j = 1:(i-1)

      k = 1
      while length(pathways[i]) > k && length(pathways[j]) > k && pathways[i][end - k] == pathways[j][end - k]
        k += 1
      end

      if debug
        println("infection of individual $(infected[i]) observed at $(obs.infectious[infected[i]])")
        println("infection pathway of individual $(infected[i]) is $(pathways[i])")
        println("infection of individual $(infected[j]) observed at $(obs.infectious[infected[j]])")
        println("infection pathway of individual $(infected[j]) is $(pathways[j])")
        if k == length(pathways[i]) || k == length(pathways[j])
          println("linear infection pathway between individual $(infected[i]) and individual $(infected[j])")
        else
          println("most recent common infection source of individuals $(infected[i]) and $(infected[j]) is $(pathways[i][end - k + 1])")
          println("the infection pathway of $(infected[i]) and $(infected[j]) diverged with $(pathways[i][end - k]) and $(pathways[j][end - k])")
        end
        println("")
      end

      if k == length(pathways[i])
        seq_dist[infected[i],infected[j]] += obs.infectious[infected[j]] - aug.exposed[pathways[j][end - k]]
        seq_dist[infected[i],infected[j]] += abs(aug.exposed[pathways[j][end - k]] - obs.infectious[infected[i]])
      elseif k == length(pathways[j])
        seq_dist[infected[i],infected[j]] += obs.infectious[infected[i]] - aug.exposed[pathways[i][end - k]]
        seq_dist[infected[i],infected[j]] += abs(aug.exposed[pathways[i][end - k]] - obs.infectious[infected[j]])
      else
        seq_dist[infected[i],infected[j]] += obs.infectious[infected[i]] - aug.exposed[pathways[i][end - k]]
        seq_dist[infected[i],infected[j]] += obs.infectious[infected[j]] - aug.exposed[pathways[j][end - k]]
        seq_dist[infected[i],infected[j]] += abs(aug.exposed[pathways[j][end - k]] - aug.exposed[pathways[i][end - k]])
      end

    end
  end

  return seq_dist += transpose(seq_dist)
end

function seq_distances(obs::SEIR_observed, aug::SEIR_augmented, infected::Vector{Int64}, network::Array, debug=false::Bool)
  """
  For a given transmission network, find the time between the pathogen sequences between every individuals i and j
  """
  pathway = [infected[1]]
  while pathway[end] != 0
    push!(pathway, findfirst(network[:,pathway[end]])-1)
  end
  pathways = Vector[pathway]

  for i = 2:length(infected)
    pathway = [infected[i]]
    while pathway[end] != 0
      push!(pathway, findfirst(network[:,pathway[end]])-1)
    end
    push!(pathways, pathway)
  end

  seq_dist = fill(0., (size(network, 2), size(network, 2)))

  for i = 1:length(infected)
    for j = 1:(i-1)

      k = 1
      while length(pathways[i]) > k && length(pathways[j]) > k && pathways[i][end - k] == pathways[j][end - k]
        k += 1
      end

      if debug
        println("infection of individual $(infected[i]) observed at $(obs.infectious[infected[i]])")
        println("infection pathway of individual $(infected[i]) is $(pathways[i])")
        println("infection of individual $(infected[j]) observed at $(obs.infectious[infected[j]])")
        println("infection pathway of individual $(infected[j]) is $(pathways[j])")
        if k == length(pathways[i]) || k == length(pathways[j])
          println("linear infection pathway between individual $(infected[i]) and individual $(infected[j])")
        else
          println("most recent common infection source of individuals $(infected[i]) and $(infected[j]) is $(pathways[i][end - k + 1])")
          println("the infection pathway of $(infected[i]) and $(infected[j]) diverged with $(pathways[i][end - k]) and $(pathways[j][end - k])")
        end
        println("")
      end

      if k == length(pathways[i])
        seq_dist[infected[i],infected[j]] += obs.infectious[infected[j]] - aug.exposed[pathways[j][end - k]]
        seq_dist[infected[i],infected[j]] += abs(aug.exposed[pathways[j][end - k]] - obs.infectious[infected[i]])
      elseif k == length(pathways[j])
        seq_dist[infected[i],infected[j]] += obs.infectious[infected[i]] - aug.exposed[pathways[i][end - k]]
        seq_dist[infected[i],infected[j]] += abs(aug.exposed[pathways[i][end - k]] - obs.infectious[infected[j]])
      else
        seq_dist[infected[i],infected[j]] += obs.infectious[infected[i]] - aug.exposed[pathways[i][end - k]]
        seq_dist[infected[i],infected[j]] += obs.infectious[infected[j]] - aug.exposed[pathways[j][end - k]]
        seq_dist[infected[i],infected[j]] += abs(aug.exposed[pathways[j][end - k]] - aug.exposed[pathways[i][end - k]])
      end

    end
  end

  return seq_dist += transpose(seq_dist)
end

# function seq_loglikelihood(seq1::Nucleotide2bitSeq, seq2::Nucleotide2bitSeq, branchdistance::Float64, substitution_matrix::Array)
#   """
#   Loglikelihood for any two aligned sequences, a specified time apart on a transmission network
#   """
#   @assert(length(seq1) == length(seq2), "Sequences not aligned")
#   ll = 0.
#   for i = 1:length(seq1)
#     base1 = convert(Int64, seq1[i])
#     for base2 = 1:4
#       if base2 == convert(Int64, seq2[i])
#         if base1 != base2
#           ll += log(1 - exp(substitution_matrix[base1, base2] .* branchdistance))
#         end
#       else
#         ll += substitution_matrix[base1, base2] .* branchdistance
#       end
#     end
#   end
#   return ll
# end

function seq_loglikelihood(seq1::Nucleotide2bitSeq, seq2::Nucleotide2bitSeq, seq_distance::Float64, substitution_matrix::Array)
  """
  Loglikelihood for any two aligned sequences, a specified time apart on a transmission network
  """
  @assert(length(seq1) == length(seq2), "Sequences not aligned")
  ll = 0.
  for i = 1:length(seq1)
    if seq1[i] == seq2[i]
      ll += -seq_distance*sum(substitution_matrix[convert(Int64, seq1[i]),:])
    else
      ll += log(1-(exp(-seq_distance*sum(substitution_matrix[convert(Int64, seq1[i]),:]))))
      ll += log(substitution_matrix[convert(Int64, seq1[i]),convert(Int64, seq2[i])]/sum(substitution_matrix[convert(Int64, seq1[i]),:]))
    end
  end
  return ll
end

function network_loglikelihood(obs::SEIR_observed, aug::SEIR_augmented, network::Array, substitution_matrix::Array, debug=false::Bool)
  """
  Loglikelihood for an entire transmission network
  """
  ll = 0.
  infected = find(isseq(obs.seq))
  seq_dist = seq_distances(obs, aug, infected, network, debug)
  for i = 1:length(infected)
    for j = 1:(i-1)
       ll += seq_loglikelihood(obs.seq[infected[i]], obs.seq[infected[j]], seq_dist[i,j], substitution_matrix)
    end
  end
  return ll
end

function SEIR_loglikelihood(α::Float64, β::Float64, η::Float64, ρ::Float64, γ::Float64, aug::SEIR_augmented, obs::SEIR_observed, dist::Metric, debug=false::Bool)
  """
  Calculate the loglikelihood and return an exposure network array under specified parameters values and observations

  α, β: powerlaw exposure kernel parameters
  η: external pressure rate
  ρ: infectivity rate (1/mean latent period)
  γ: recovery rate (1/mean infectious period)
  ν: detection rate (1/mean detection lag)
  """
  # Initiate an exposure network
  network = fill(false, (1 + length(obs.covariates), length(obs.covariates)))

  ll = loglikelihood(Exponential(1/γ), (aug.removed .- aug.infectious)[!isnan(aug.removed)])

  if debug && ll == -Inf
    print("Infectious period caused loglikelihood to go to -Inf")
  end

  # Create event timing array
  event_times = [aug.exposed aug.infectious aug.removed]

  # Find event order
  event_order = sortperm(event_times[:])

  # Create empty rate array
  rate_array = fill(0., (1 + length(obs.covariates) + 2, length(obs.covariates)))

  # First row is external pressure rate
  rate_array[1,:] = η

  # The rest of the events
  for i = 1:length(event_order)
    isnan(event_times[event_order[i]]) && break
    rate_array_sum = sum(rate_array)
    id = ind2sub(size(event_times), event_order[i])

    # Don't consider likelilihood contribution of first event
    if i > 1
      ll += logpdf(Exponential(1/rate_array_sum), event_times[event_order[i]] - event_times[event_order[i-1]])
      ll += log(sum(rate_array[:,id[1]])/rate_array_sum)
    end

    # Exposure event
    if id[2] == 1
      # Generate a exposure source based on disease pressures at time of exposure
      network[:, id[1]] = rand(Multinomial(1, rate_array[1:(length(obs.covariates)+1), id[1]]./sum(rate_array[1:(length(obs.covariates)+1), id[1]])))
      # Update rate array (exposure rates, latent period)
      rate_array[:, id[1]] = 0.
      rate_array[1 + size(rate_array, 2) + 1, id[1]] = ρ

    # Infectiousness event
    elseif id[2] == 2
      # Update rate_array for latent and infectious periods
      rate_array[1 + size(rate_array, 2) + 1, id[1]] = 0.
      rate_array[1 + size(rate_array, 2) + 2, id[1]] = γ
      # Update exposure rates for rest of susceptible population
      for j = 1:size(rate_array, 2)
        if j != id[1] && rate_array[1, j] != 0.
          rate_array[id[1] + 1, j] = α*evaluate(dist, obs.covariates[id[1]], obs.covariates[j])^-β
#         else
#           rate_array[id[1] + 1, j] = 0.
        end
      end

    # Removal event
    elseif id[2] == 3
      # Update rate_array for exposure & removal
      rate_array[id[1] + 1, :] = 0.
      rate_array[1 + size(rate_array,2) + 2, id[1]] = 0.
    end

    # Provide loop position when loglikelihood goes to -Inf if desired
    if debug && ll == -Inf
      if id[2] == 1
        print("Event $i (exposure of individual $(id[1])) caused loglikelihood to go to -Inf")
      elseif id[2] == 2
        print("Event $i (infection of individual $(id[1])) caused loglikelihood to go to -Inf")
      elseif id[2] == 3
        print("Event $i (removal of individual $(id[1])) caused loglikelihood to go to -Inf")
      end
    end

    # Prevent needless calculation
    if ll == -Inf
      break
    end

  end
  return ll, network
end

function initialize(ilm_priors::SEIR_priors, mutation_priors::JC69_priors, detection_priors::Lag_priors, obs::SEIR_observed, limit=500::Int, debug=false::Bool, dist=Euclidean())
  """
  Initiate an Trace object by sampling from specified prior distributions
  """
  ilm_params = randprior(ilm_priors)
  mutation_params = randprior(mutation_priors)
  detection_params = randprior(detection_priors)
  aug = augment(ilm_params[4], detection_params[1], obs)
  ll, network = SEIR_loglikelihood(ilm_params[1], ilm_params[2], ilm_params[3], ilm_params[4], ilm_params[5], aug, obs, dist, debug)
  ll += network_loglikelihood(obs, aug, network, jc69([mutation_params[1]]), debug)
  count = 1

  # Retry initialization until non-negative infinity loglikelihood
  while ll == -Inf && count < limit
    count += 1
    ilm_params = randprior(ilm_priors)
    mutation_params = randprior(mutation_priors)
    detection_params = randprior(detection_priors)
    aug = augment(ilm_params[4], detection_params[1], obs)
    ll, network = SEIR_loglikelihood(ilm_params[1], ilm_params[2], ilm_params[3], ilm_params[4], ilm_params[5], aug, obs, dist, debug)
    ll += network_loglikelihood(obs, aug, network, jc69([mutation_params[1]]), debug)
  end

  if count < limit
    print("Successfully initialized on attempt $count")
    lp = ll + logprior(ilm_priors, ilm_params) + logprior(mutation_priors, mutation_params) + logprior(detection_priors, detection_params)
      return SEIR_trace([ilm_params[1]], [ilm_params[2]], [ilm_params[3]], [ilm_params[4]], [ilm_params[5]], [aug], Array[network], [lp]), Lag_trace([detection_params[1]]), JC69_trace([mutation_params[1]])
  else
    print("Failed to initialize after $count attempts")
  end
end

function initialize(ilm_priors::SEIR_priors, mutation_priors::JC69_priors, obs::SEIR_observed, limit=500::Int, debug=false::Bool, dist=Euclidean())
  """
  Initiate an Trace object by sampling from specified prior distributions
  """
  ilm_params = randprior(ilm_priors)
  mutation_params = randprior(mutation_priors)
  aug = augment(ilm_params[4], obs)
  ll, network = SEIR_loglikelihood(ilm_params[1], ilm_params[2], ilm_params[3], ilm_params[4], ilm_params[5], aug, obs, dist, debug)
  ll += network_loglikelihood(obs, aug, network, jc69([mutation_params[1]]), debug)
  count = 1

  # Retry initialization until non-negative infinity loglikelihood
  while ll == -Inf && count < limit
    count += 1
    ilm_params = randprior(ilm_priors)
    mutation_params = randprior(mutation_priors)
    aug = augment(ilm_params[4], obs)
    ll, network = SEIR_loglikelihood(ilm_params[1], ilm_params[2], ilm_params[3], ilm_params[4], ilm_params[5], aug, obs, dist, debug)
    ll += network_loglikelihood(obs, aug, network, jc69([mutation_params[1]]), debug)
  end

  if count < limit
    print("Successfully initialized on attempt $count")
    lp = ll + logprior(ilm_priors, ilm_params) + logprior(mutation_priors, mutation_params)
      return SEIR_trace([ilm_params[1]], [ilm_params[2]], [ilm_params[3]], [ilm_params[4]], [ilm_params[5]], [aug], Array[network], [lp]), JC69_trace([mutation_params[1]])
  else
    print("Failed to initialize after $count attempts")
  end
end

function initialize(ilm_priors::SEIR_priors, detection_priors::Lag_priors, obs::SEIR_observed, limit=500::Int, debug=false::Bool, dist=Euclidean())
  """
  Initiate an Trace object by sampling from specified prior distributions
  """
  ilm_params = randprior(ilm_priors)
  detection_params = randprior(detection_priors)
  aug = augment(ilm_params[4], detection_params[1], obs)
  ll, network = SEIR_loglikelihood(ilm_params[1], ilm_params[2], ilm_params[3], ilm_params[4], ilm_params[5], aug, obs, dist, debug)
  count = 1

  # Retry initialization until non-negative infinity loglikelihood
  while ll == -Inf && count < limit
    count += 1
    ilm_params = randprior(ilm_priors)
    detection_params = randprior(detection_priors)
    aug = augment(ilm_params[4], detection_params[1], obs)
    ll, network = SEIR_loglikelihood(ilm_params[1], ilm_params[2], ilm_params[3], ilm_params[4], ilm_params[5], aug, obs, dist, debug)
  end

  if count < limit
    print("Successfully initialized on attempt $count")
    lp = ll + logprior(ilm_priors, ilm_params) + logprior(detection_priors, detection_params)
      return SEIR_trace([ilm_params[1]], [ilm_params[2]], [ilm_params[3]], [ilm_params[4]], [ilm_params[5]], [aug], Array[network], [lp]), Lag_trace([detection_params[1]])
  else
    print("Failed to initialize after $count attempts")
  end
end

function initialize(ilm_priors::SEIR_priors, obs::SEIR_observed, limit=500::Int, debug=false::Bool, dist=Euclidean())
  """
  Initiate an Trace object by sampling from specified prior distributions
  """
  ilm_params = randprior(ilm_priors)
  aug = augment(ilm_params[4], obs)
  ll, network = SEIR_loglikelihood(ilm_params[1], ilm_params[2], ilm_params[3], ilm_params[4], ilm_params[5], aug, obs, dist, debug)
  count = 1

  # Retry initialization until non-negative infinity loglikelihood
  while ll == -Inf && count < limit
    count += 1
    ilm_params = randprior(ilm_priors)
    aug = augment(ilm_params[4], obs)
    ll, network = SEIR_loglikelihood(ilm_params[1], ilm_params[2], ilm_params[3], ilm_params[4], ilm_params[5], aug, obs, dist, debug)
  end

  if count < limit
    print("Successfully initialized on attempt $count")
    lp = ll + logprior(ilm_priors, ilm_params)
      return SEIR_trace([ilm_params[1]], [ilm_params[2]], [ilm_params[3]], [ilm_params[4]], [ilm_params[5]], [aug], Array[network], [lp])
  else
    print("Failed to initialize after $count attempts")
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

  @assert(size(transition_cov) == (7,7), "transition_cov must be a 7x7 matrix")

  for i = 1:n

    # Create and incremenet progress bar
    if progress
      if i == 1
        progressbar = Progress(n, 5, "Performing $n MCMC iterations...", 30)
      else
        next!(progressbar)
      end
    end

    # Only generate valid proposals
    step = rand(MvNormal(transition_cov))
    ilm_proposal = [ilm_trace.α[end], ilm_trace.β[end], ilm_trace.η[end], ilm_trace.ρ[end], ilm_trace.γ[end]] .+ step[1:5]
    detection_proposal = [detection_trace.ν[end]] .+ step[6]
    mutation_proposal = [mutation_trace.λ[end]] .+ step[7]

    lp = logprior(ilm_priors, ilm_proposal)
    lp += logprior(detection_priors, detection_proposal)
    lp += logprior(mutation_priors, mutation_proposal)

    while lp == -Inf
      step = rand(MvNormal(transition_cov))
      ilm_proposal = [ilm_trace.α[end], ilm_trace.β[end], ilm_trace.η[end], ilm_trace.ρ[end], ilm_trace.γ[end]] .+ step[1:5]
      detection_proposal = [detection_trace.ν[end]] .+ step[6]
      mutation_proposal = [mutation_trace.λ[end]] .+ step[7]

      lp = logprior(ilm_priors, ilm_proposal)
      lp += logprior(detection_priors, detection_proposal)
      lp += logprior(mutation_priors, mutation_proposal)
    end

    # Augment the data
    aug = augment(ilm_proposal[4], detection_proposal[1], obs)

    # ILM loglikelihood component
    ll, network = SEIR_loglikelihood(ilm_proposal[1], ilm_proposal[2], ilm_proposal[3], ilm_proposal[4], ilm_proposal[5], aug, obs, dist)
    lp += ll

    # Network loglikelihood
    if lp > -Inf
      lp += network_loglikelihood(obs, aug, network, jc69([mutation_proposal[1]]), debug)
    end

    # Accept/reject based on logposterior
    if lp == -Inf
      reject = true
    elseif lp > ilm_trace.logposterior[end]
      reject = false
    elseif exp(lp - ilm_trace.logposterior[end]) >= rand()
      reject = false
    else
      reject = true
    end

    if reject
      ilm_proposal[1] = ilm_trace.α[end]
      ilm_proposal[2] = ilm_trace.β[end]
      ilm_proposal[3] = ilm_trace.η[end]
      ilm_proposal[4] = ilm_trace.ρ[end]
      ilm_proposal[5] = ilm_trace.γ[end]
      aug = ilm_trace.aug[end]
      network = ilm_trace.network[end]
      lp = ilm_trace.logposterior[end]

      detection_proposal[1] = detection_trace.ν[end]

      mutation_proposal[1] = mutation_trace.λ[end]

#       # Re-augment data (and recalculate log posterior) if proposal is rejected...
#       aug = augment(ilm_proposal[4], detection_proposal[1], obs)
#       lp, network = SEIR_loglikelihood(ilm_proposal[1], ilm_proposal[2], ilm_proposal[3], ilm_proposal[4], ilm_proposal[5], aug, obs, dist)
#       lp += network_loglikelihood(obs, aug, network, jc69([mutation_proposal[1]]), debug)
#       lp += logprior(ilm_priors, ilm_proposal)
#       lp += logprior(detection_priors, detection_proposal)
#       lp += logprior(mutation_priors, mutation_proposal)
    end

    # Update trace objects
    push!(ilm_trace.α, ilm_proposal[1])
    push!(ilm_trace.β, ilm_proposal[2])
    push!(ilm_trace.η, ilm_proposal[3])
    push!(ilm_trace.ρ, ilm_proposal[4])
    push!(ilm_trace.γ, ilm_proposal[5])
    push!(ilm_trace.aug, aug)
    push!(ilm_trace.network, network)
    push!(ilm_trace.logposterior, lp)

    push!(detection_trace.ν, detection_proposal[1])

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

  @assert(size(transition_cov) == (7,7), "transition_cov must be a 7x7 matrix")

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
