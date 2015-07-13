"""
infer.jl - pathogen evolution and transmission dynamic inference tools
Justin Angevaare
June 2015
"""

function SEIR_surveilance(population::Population, ν::Float64)
  """
  Gather surveillance data on specific individuals in a population, with an exponentially distributed detection lag with rate ν
  """
  @assert(0. < ν, "ν, the detection rate parameter must be greater than 0")
  @warning("Movement and covariate changes currently not supported, only initial conditions considered")
  exposed_actual = fill(NaN, length(ids))
  infectious_actual = fill(NaN, length(ids))
  infectious_observed = fill(NaN, length(ids))
  removed_actual = fill(NaN, length(ids))
  removed_observed = fill(NaN, length(ids))
  covariates = fill(fill(NaN, length(population.history[ids[1]][1][1])), length(ids))
  seq = fill(NaN, length(ids))

  for i = 2:length(population.events)
    # Initial conditions
    covariates[i] = population.history[i][1][1]

    # Exposure time (unobservable)
    if length(population.events[i][1]) > 0
      exposed_actual[i] = population.events[i][1][1]
    end

    # Infectious time (observed with latency)
    if length(population.events[i][3]) > 0
      infectious_actual[i] = population.events[i][3][1]
      infectious_observed[i] = infectious_actual[i] + rand(Exponential(1/ν))
      seq[i] = population.history[i][2][find(infectious_observed[i] .<= population.events[i][6])[end]]
    end

    # Removal time (observed with latency)
    if length(population.events[i][4]) > 0
      removed_actual[i] = population.events[i][4][1]
      removed_observed[i] = removed_actual[i] + rand(Exponential(1/ν))
    end
  end
  return SEIR_events(exposed, infectious_actual, infectious_observed, removed_actual, removed_observed, covariates, seq)
end

function SEIR_augmentation(ρ::Float64, ν::Float64, obs::SEIR_events)
  """
  Augments surveilance data, organizes observations
  """
  infectious_augmented = SEIR_events.infectious_observed - rand(Exponential(1/ν), length(SEIR_events.infectious_observed))
  exposed_augmented = infectious_augmented - rand(Exponential(1/ρ), length(SEIR_events.infectious_observed))
  recovered_augmented = SEIR_events.recovered_observed - rand(Exponential(1/ν), length(SEIR_events.recovered_observed))
  return SEIR_augmented(infectious_augmented, exposed_augmented, recovered_augmented)
end

function SEIR_logprior(α_prior::UnivariateDistribution, β_prior::UnivariateDistribution, ρ_prior::UnivariateDistribution, γ_prior::UnivariateDistribution, η_prior::UnivariateDistribution, ν_prior::UnivariateDistribution)
  """
  Create a log prior function using specified prior univariate distributions
  α, β: powerlaw exposure kernel parameters
  η: external pressure rate
  ρ: infectivity rate (1/mean latent period)
  γ: recovery rate (1/mean infectious period)
  ν: detection rate (1/mean detection lag)
  """
  return function(α::Float64, β::Float64, ρ::Float64, γ::Float64, η::Float64, ν::Float64)
    """
    α, β: powerlaw exposure kernel parameters
    η: external pressure rate
    ρ: infectivity rate (1/mean latent period)
    γ: recovery rate (1/mean infectious period)
    ν: detection rate (1/mean detection lag)
    """
    return logpdf(α_prior, α) + logpdf(β_prior, β) + logpdf(ρ_prior, ρ) + logpdf(γ_prior, γ) + logpdf(η_prior, η) + logpdf(ν_prior, ν)
  end
end

function SEIR_loglikelihood(α::Float64, β::Float64, ρ::Float64, γ::Float64, η::Float64, ν::Float64, aug::SEIR_augmented, obs::SEIR_events, dist=Euclidean())
  """
  α, β: powerlaw exposure kernel parameters
  η: external pressure rate
  ρ: infectivity rate (1/mean latent period)
  γ: recovery rate (1/mean infectious period)
  ν: detection rate (1/mean detection lag)
  """
  ll = loglikelihood(Exponential(1/γ), (aug.recovered_augmented .- aug.infectious_augmented)[!isnan(aug.recovered_augmented)])
# DO NOT CONSIDERED LIKELIHOOD OF PURELY AUGMENTED DATA
#  ll = loglikelihood(Exponential(1/ρ), aug.infectious_augmented .- aug.exposed_augmented)
#  ll += loglikelihood(Exponential(1/ν), obs.recovered_observed .- aug.recovered_augmented)
#  ll += loglikelihood(Exponential(1/ν), obs.infectious_observed .- aug.infectious_augmented)

  # Create event timing array
  event_times = [aug.exposed_augmented
                 aug.infectious_augmented
                 recovered_augmented]

  # Find event order
  event_order = sortperm(event_times[:])

  # Create empty rate array
  rate_array = fill(0., (1 + length(obs.id) + 2, length(obs.id)))

  # First row is external pressure rate
  rate_array[1,] = η

  # The rest of the events
  for i = 1:length(event_order)
    isnan(event_times[event_order[i]]) && break

    id = ind2sub(event_order[i], size(event_times))
    ll += loglikelihood(Exponential(1/rate_array_sum), event_times[event_order[i]])
    ll += log(sum(rate_array[:,id[1]])/sum(rate_array))

    # Update rate_array for latent period
    if id[2] == 1
      rate_array[:, id[1]] = 0.
      rate_array[size(rate_array,2) + 1, id[1]] = ρ
    end

    # Update rate_array for infectious period & exposure
    if id[2] == 2
      rate_array[size(rate_array, 2) + 1, id[1]] = 0.
      rate_array[size(rate_array, 2) + 2, id[1]] = γ
      for j = size(rate_array, 2)
        if j != id[1]
          rate_array[id[1], j] = α*evaluate(dist, obs.covariates[id[1]], obs.covariates[j])^-β
        else
          rate_array[id[1], j] = 0.
        end
      end
    end

    # Update rate_array for recovery & exposure
    if id[2] == 3
      rate_array[id[1],:] = 0.
      rate_array[size(rate_array,2) + 2, id[1]] = 0.
    end
  end
  return ll
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
