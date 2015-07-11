"""
infer.jl - pathogen evolution and transmission dynamic inference tools
Justin Angevaare
June 2015
"""

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

function treedistance(leaf1::Int64, leaf2::Int64, tree::Tree)
  """
  Find the minimumum branch distance between two leaves
  """
  @assert(all(1 .<= [leaf1, leaf2] .<= length(tree.distances)), "Invalid leaves specified")
  depthlimit = minimum([length(tree.positions[leaf1]), length(tree.positions[leaf2])])
  sharednode = findfirst(tree.positions[leaf1][1:depthlimit] .!= tree.positions[leaf2][1:depthlimit])
  return sum([tree.distances[leaf1][sharednode:end], tree.distances[leaf2][sharednode:end]])
end

function seqdistance(seq1::Nucleotide2bitSeq, seq2::Nucleotide2bitSeq, substitution_matrix::Array)
  """
  Find the distance between two sequences as per a specificied substitution rate matrix
  """
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

function augment1(ρ::Float64, ν::Float64, obs)
  """
  Augments surveilance data with infectiousness
  """
  returns aug
  end

function loglikelihood1(α::Float64, β::Float64, ρ::Float64, γ::Float64, η::Float64, ν::Float64, obs, aug)
  """
  α, β: powerlaw exposure kernel parameters
  η: external pressure rate
  ρ: infectivity rate (1/mean latent period)
  γ: recovery rate (1/mean infectious period)
  ν: detection rate (1/mean detection lag)
  """
end

