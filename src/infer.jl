"""
infer.jl - pathogen evolution and transmission dynamic inference tools
Justin Angevaare
June 2015
"""

function SEIR_surveilance(ids::Vector{Int64}, population::Population, ν::Float64)
  """
  Gather surveillance data on specific individuals in a population, with an exponentially distributed detection lag with rate ν
  """
  @assert(all(2 .<= ids .<= length(population.events)), "Invalid ID provided")
  @assert(0. < ν, "ν, the detection rate parameter must be greater than 0")
  @warning("Movement and covariate changes currently not supported, only initial conditions considered")
  exposed = fill(NaN, length(ids))
  infectious_actual = fill(NaN, length(ids))
  infectious_observed = fill(NaN, length(ids))
  removed_actual = fill(NaN, length(ids))
  removed_observed = fill(NaN, length(ids))
  covariates = fill(fill(NaN, length(population.history[ids[1]][1][1])), length(ids))
  seq = fill(NaN, length(ids))

  for i = 1:length(ids)
    # Initial conditions
    covariates[i] = population.history[ids[i]][1][1]

    # Exposure time (unobservable)
    if length(population.events[ids[i]][1]) > 0
      exposed[i] = population.events[ids[i]][1][1]
    end

    # Infectious time (observed with latency)
    if length(population.events[ids[i]][3]) > 0
      infectious_actual[i] = population.events[ids[i]][3][1]
      infectious_observed[i] = infectious_actual[i] + rand(Exponential(1/ν))
      seq[i] = population.history[ids[i]][2][find(infectious_observed[i] .<= population.events[ids[i]][6])[end]]
    end

    # Removal time (observed with latency)
    if length(population.events[ids[i]][4]) > 0
      removed_actual[i] = population.events[ids[i]][4][1]
      removed_observed[i] = removed_actual[i] + rand(Exponential(1/ν))
    end
  end
  return SEIR_events(exposed, infectious_actual, infectious_observed, removed_actual, removed_observed, covariates, seq)
end

function create_tree(sequences::Vector{Nucleotide2bitSeq}, times::Vector{Float64})
  """
  Generate a phylogenetic tree based on sample times and sequences
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

function create_logprior1(α_prior::UnivariateDistribution, β_prior::UnivariateDistribution, ρ_prior::UnivariateDistribution, γ_prior::UnivariateDistribution, η_prior::UnivariateDistribution, ν_prior::UnivariateDistribution)
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

function augorg1(ρ::Float64, ν::Float64, obs::DataFrame)
  """
  Augments surveilance data, organizes observations
  """
  augorg = [obs, DataFrame(detectionlag = rand(Exponential(1/ν), size(obs, 1)))]
  augorg[:truetime] = augorg[:time] - augorg[:detectionlag]

  return augorg
  end

function loglikelihood1(α::Float64, β::Float64, ρ::Float64, γ::Float64, η::Float64, ν::Float64, augorg::DataFrame)
  """
  α, β: powerlaw exposure kernel parameters
  η: external pressure rate
  ρ: infectivity rate (1/mean latent period)
  γ: recovery rate (1/mean infectious period)
  ν: detection rate (1/mean detection lag)
  """
  loglikelihood(Exponential(1/ρ), aug.latentperiod[end]) + loglikelihood(Exponential(1/γ), aug.infectiousperiod[end]) + loglikelihood(Exponential(1/ν), aug.detectionlag[end])
  return

end

