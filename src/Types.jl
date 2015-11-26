import Base.display
import Base.push!

"""
A specialized type describing individuals in a population
"""
type Population
  events::Array
  history::Array
  timeline::Array
end


"""
Contains rates, and information to as which event they refer to
"""
type RateArray
  rates::Array
  events::Array
end

abstract Actual

abstract Observed

abstract Augmented

"""
Contains the actual event times
"""
type SEIR_actual <: Actual
  exposed::Vector{Float64}
  infectious::Vector{Float64}
  removed::Vector{Float64}
  covariates::Vector{Vector{Float64}}
  seq::Vector{Any}
end


"""
Contains event observations
"""
type SEIR_observed <: Observed
  infectious::Vector{Float64}
  removed::Vector{Float64}
  covariates::Vector{Vector{Float64}}
  seq::Vector{Any}
  sequenced::Vector{Bool}
end


"""
Contains augmented event times
"""
type SEIR_augmented <: Augmented
  exposed::Vector{Float64}
  infectious::Vector{Float64}
  removed::Vector{Float64}
end


"""
Contains the actual event times
"""
type SIR_actual <: Actual
  infectious::Vector{Float64}
  removed::Vector{Float64}
  covariates::Vector{Vector{Float64}}
  seq::Vector{Any}
end


"""
Contains event observations
"""
type SIR_observed <: Observed
  infectious::Vector{Float64}
  removed::Vector{Float64}
  covariates::Vector{Vector{Float64}}
  seq::Vector{Any}
  sequenced::Vector{Bool}
end


"""
Contains augmented event times
"""
type SIR_augmented <: Augmented
  infectious::Vector{Float64}
  removed::Vector{Float64}
end


abstract Priors


abstract ILM_priors <: Priors


abstract Detection_priors <: Priors


abstract Mutation_priors <: Priors


"""
Prior distributions for SEIR model inference
"""
type SEIR_priors <: ILM_priors
  α::UnivariateDistribution
  β::UnivariateDistribution
  η::UnivariateDistribution
  ρ::UnivariateDistribution
  γ::UnivariateDistribution
end


"""
Prior distributions for SEIR model inference
"""
type SIR_priors <: ILM_priors
  α::UnivariateDistribution
  β::UnivariateDistribution
  η::UnivariateDistribution
  γ::UnivariateDistribution
end


"""
Prior distributions for detection lag model inference
"""
type Lag_priors{T<:UnivariateDistribution} <: Detection_priors
  ν::T
end


"""
Prior distributions for JC69 model inference
"""
type JC69_priors{T<:UnivariateDistribution} <: Mutation_priors
  λ::T
end


abstract Trace


abstract ILM_trace <: Trace


abstract Detection_trace <: Trace


abstract Mutation_trace <: Trace


"""
Contains an MCMC trace for an SEIR model
"""
type SEIR_trace <: ILM_trace
  α::Vector{Float64}
  β::Vector{Float64}
  η::Vector{Float64}
  ρ::Vector{Float64}
  γ::Vector{Float64}
  aug::Vector{SEIR_augmented}
  network_rates::Vector{Array{Float64}}
  network::Vector{Array{Bool}}
  logposterior::Vector{Float64}
end


"""
Contains an MCMC trace for an SIR model
"""
type SIR_trace <: ILM_trace
  α::Vector{Float64}
  β::Vector{Float64}
  η::Vector{Float64}
  γ::Vector{Float64}
  aug::Vector{SIR_augmented}
  network_rates::Vector{Array{Float64}}
  network::Vector{Array{Bool}}
  logposterior::Vector{Float64}
end


"""
Contains an MCMC trace object for detection lag model
"""
type Lag_trace <: Detection_trace
  ν::Vector{Float64}
end


"""
Contains an MCMC trace object for JC69 model parameters
"""
type JC69_trace <: Mutation_trace
  λ::Vector{Float64}
end
