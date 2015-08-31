"""
types.jl
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

abstract Priors
abstract Mutation_priors <: Priors

type ILM_priors{T<:UnivariateDistribution} <: Priors
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

type JC69_priors{T<:UnivariateDistribution} <: Mutation_priors
  """
  Prior distributions for JC69 model inference
  """
  λ::T
end

abstract Trace
abstract Mutation_trace <: Trace

type ILM_trace <: Trace
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

type JC69_trace <: Mutation_trace
  """
  Contains an MCMC trace object for JC69 model parameters
  """
  λ::Vector{Float64}
end
