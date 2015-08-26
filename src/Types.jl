"""
types.jl
Justin Angevaare
"""

import Base.display
import Base.push!

type Population
  """
  A specialized type where each individual in a population is characterized by a vector of vectors for event details and by a vector of arrays for event histories
  """
  events::Array
  history::Array
  timeline::Array
end

type RateArray
  """
  Contains rates, and information to as what they refer to
  """
  rates::Array
  events::Array
end

type SEIR_actual
  """
  Contains actual event times
  """
  exposed::Vector{Float64}
  infectious::Vector{Float64}
  removed::Vector{Float64}
  covariates::Vector{Vector{Float64}}
  seq::Vector{Any}
end

type SEIR_observed
  """
  Contains observed event times and information
  """
  infectious::Vector{Float64}
  removed::Vector{Float64}
  covariates::Vector{Vector{Float64}}
  seq::Vector{Any}
end

type SEIR_augmented
  """
  Contains event times from data augmentation
  """
  exposed::Vector{Float64}
  infectious::Vector{Float64}
  removed::Vector{Float64}
end

type SEIR_priors{T<:UnivariateDistribution}
  """
  Prior distributions for SEIR model inference
  """
  α::T
  β::T
  η::T
  ρ::T
  γ::T
  ν::T
end

type PhyloSEIR_priors{T<:UnivariateDistribution}
  """
  Prior distributions for SEIR model inference
  """
  α::T
  β::T
  η::T
  ρ::T
  γ::T
  ν::T
  mutation::Tuple{T}
end

type SEIR_trace
  """
  Contains an MCMC trace object
  """
  α::Vector{Float64}
  β::Vector{Float64}
  η::Vector{Float64}
  ρ::Vector{Float64}
  γ::Vector{Float64}
  ν::Vector{Float64}
  aug::Vector{SEIR_augmented}
  network::Vector{Array{Bool}}
  logposterior::Vector{Float64}
end

type PhyloSEIR_trace
  """
  Contains an MCMC trace object
  """
  α::Vector{Float64}
  β::Vector{Float64}
  η::Vector{Float64}
  ρ::Vector{Float64}
  γ::Vector{Float64}
  ν::Vector{Float64}
  mutation::Array{Float64}
  aug::Vector{SEIR_augmented}
  network::Vector{Array{Bool}}
  logposterior::Vector{Float64}
end
