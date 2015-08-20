"""
infer.jl
Pathogen evolution and transmission dynamic inference tools
Justin Angevaare
"""

function SEIR_surveilance(population::Population, ν::Float64)
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

function SEIR_augmentation(ρ::Float64, ν::Float64, obs::SEIR_observed)
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

function SEIR_loglikelihood(α::Float64, β::Float64, η::Float64, ρ::Float64, γ::Float64, ν::Float64, aug::SEIR_augmented, obs::SEIR_observed, dist::Metric, debug=false::Bool)
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

function SEIR_logprior(priors::SEIR_priors, α::Float64, β::Float64, η::Float64, ρ::Float64, γ::Float64, ν::Float64)
  """
  Calculate the logprior from prior distributions defined in `SEIR_priors` and specific parameter values
  """
  if ν < Inf
    return logpdf(priors.α, α) + logpdf(priors.β, β) + logpdf(priors.η, η) + logpdf(priors.ρ, ρ) + logpdf(priors.γ, γ)  + logpdf(priors.ν, ν)
  elseif ν == Inf
    return logpdf(priors.α, α) + logpdf(priors.β, β) + logpdf(priors.η, η) + logpdf(priors.ρ, ρ) + logpdf(priors.γ, γ)
  end
end

function SEIR_initialize(priors::SEIR_priors, obs::SEIR_observed, limit=1000::Int, debug=false::Bool, dist=Euclidean())
  """
  Initiate an SEIR_trace object by sampling from specified prior distributions

  α, β: powerlaw exposure kernel parameters
  η: external pressure rate
  ρ: infectivity rate (1/mean latent period)
  γ: recovery rate (1/mean infectious period)
  ν: detection rate (1/mean detection lag)
  """
  α = rand(priors.α)
  β = rand(priors.β)
  η = rand(priors.η)
  ρ = rand(priors.ρ)
  γ = rand(priors.γ)
  ν = rand(priors.ν)
  aug = SEIR_augmentation(ρ, ν, obs)
  ll, network = SEIR_loglikelihood(α, β, η, ρ, γ, ν, aug, obs, dist, debug)
  count = 1

  # Retry initialization until non-negative infinity loglikelihood
  while ll == -Inf && count < limit
    count += 1
    α = rand(priors.α)
    β = rand(priors.β)
    η = rand(priors.η)
    ρ = rand(priors.ρ)
    γ = rand(priors.γ)
    ν = rand(priors.ν)
    aug = SEIR_augmentation(ρ, ν, obs)
    ll, network = SEIR_loglikelihood(α, β, η, ρ, γ, ν, aug, obs, dist)
  end

  if count < limit
    print("Successfully initialized on attempt $count")
    logposterior = ll + SEIR_logprior(priors, α, β, η, ρ, γ, ν)
    return SEIR_trace([α], [β], [η], [ρ], [γ], [ν], [aug], Array[network], [logposterior])
  else
    print("Failed to initialize after $count attempts")
  end
end

function SEIR_MCMC(n::Int64, transition_cov::Array{Float64}, trace::SEIR_trace, priors::SEIR_priors, obs::SEIR_observed, dist=Euclidean())
  """
  Performs `n` data-augmented metropolis hastings MCMC iterations. Initiates a single chain by sampling from prior distribution

  α, β: powerlaw exposure kernel parameters
  η: external pressure rate
  ρ: infectivity rate (1/mean latent period)
  γ: recovery rate (1/mean infectious period)
  ν: detection rate (1/mean detection lag)
  """
  @assert(size(transition_cov) == (6,6), "transition_cov must be a 6x6 matrix")

  for i = 1:n

    # Create and incremenet progress bar
    if i == 1
      progressbar = Progress(n, 5, "Performing $n MCMC iterations...", 30)
    else
      next!(progressbar)
    end

    # Only generate valid proposals
    proposal = [trace.α[end], trace.β[end], trace.η[end], trace.ρ[end], trace.γ[end], trace.ν[end]] .+ rand(MvNormal(transition_cov))
    while any(proposal .< 0.)
      proposal = [trace.α[end], trace.β[end], trace.η[end], trace.ρ[end], trace.γ[end], trace.ν[end]] .+ rand(MvNormal(transition_cov))
    end

    # Augment the data
    aug = SEIR_augmentation(proposal[4], proposal[6], obs)

    # Loglikelihood calculation and exposure network array
    ll, network = SEIR_loglikelihood(proposal[1], proposal[2], proposal[3], proposal[4], proposal[5], proposal[6], aug, obs, dist)

    # Add logprior for logposterior
    logposterior = ll + SEIR_logprior(priors, proposal[1], proposal[2], proposal[3], proposal[4], proposal[5], proposal[6])

    # Accept/reject based on logposterior
    if logposterior == -Inf
      reject = true
    elseif logposterior > trace.logposterior[end]
      reject = false
    elseif exp(logposterior - trace.logposterior[end]) >= rand()
      reject = false
    else
      reject = true
    end

    # Re-augment data (and recalculate log posterior) if proposal is rejected...
    if reject
      proposal[1] = trace.α[end]
      proposal[2] = trace.β[end]
      proposal[3] = trace.η[end]
      proposal[4] = trace.ρ[end]
      proposal[5] = trace.γ[end]
      proposal[6] = trace.ν[end]

      aug = SEIR_augmentation(proposal[4], proposal[6], obs)
      ll, network = SEIR_loglikelihood(proposal[1], proposal[2], proposal[3], proposal[4], proposal[5], proposal[6], aug, obs, dist)
      logposterior = ll + SEIR_logprior(priors, proposal[1], proposal[2], proposal[3], proposal[4], proposal[5], proposal[6])
    end

    # Update chain
    push!(trace.α, proposal[1])
    push!(trace.β, proposal[2])
    push!(trace.η, proposal[3])
    push!(trace.ρ, proposal[4])
    push!(trace.γ, proposal[5])
    push!(trace.ν, proposal[6])
    push!(trace.aug, aug)
    push!(trace.network, network)
    push!(trace.logposterior, logposterior)
  end

  return trace
end

function create_tree(sequences::Vector{Nucleotide2bitSeq}, times::Vector{Float64})
  """
  Generate a phylogenetic tree based on sample times and sequences
  To do: update to utilize SEIR_events
  """
  @assert(length(sequences)==length(times), "There must be one sample time for each sequence")
  @assert(length(sequences)>2, "There must be at least 3 samples")
  # root
  vertices = TreeVertex()
  # nodes
  for i = 1:(length(sequences) - 2)
    push!(vertices, TreeVertex(minimum(times)))
  end
  # leaves
  for i = 1:length(sequences)
    push!(vertices, TreeVertex(sequences[i], times[i]))
  end
  # Create edges
  edges = Vector{TreeEdge}
  for i = 1:length(vertices)
    for j = 1:length(vertices)
      if vertices[i].out & vertices[j].in
        push!(edges, TreeEdge(i, j))
      end
    end
  end
  return Tree(vertices, edges)
end

function seqdistance(ancestor::Nucleotide2bitSeq, descendent::Nucleotide2bitSeq, substitution_matrix::Array)
  """
  Compute the genetic distance between two nucleotide sequences based on a `substitution_matrix`
  """
  @assert(length(ancestor) == length(descendent), "Sequences must be equal in length")
  rate_vector = Float64[]
  for i = 1:length(ancestor)
    if ancestor[i] != descendent[i]
     push!(rate_vector, substitution_matrix[convert(Int64, ancestor[i]), convert(Int64, descendent[i])])
    end
  end
  rate_vector .^= -1
  return sum(rate_vector)
end

function branchloglikelihood(seq1::Nucleotide2bitSeq, seq2::Nucleotide2bitSeq, branchdistance::Float64, substitution_matrix::Array)
  """
  Log likelihood for any two aligned sequences, a specified distance apart on a phylogenetic tree
  """
  @assert(length(seq1) == length(seq2), "Sequences not aligned")
  ll = 0
  for i = 1:length(seq1)
    base1 = convert(Int64, seq1[i])
    for base2 = 1:4
      if base2 == convert(Int64, seq2[i])
        if base1 != base2
          ll += log(1 - exp(substitution_matrix[base1, base2] .* branchdistance))
        end
      else
        ll += substitution_matrix[base1, base2] .* branchdistance
      end
    end
  end
  return ll
end
