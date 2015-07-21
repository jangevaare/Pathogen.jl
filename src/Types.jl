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

type TreeVertex
  """
  Describes a vertex in a phylogenetic tree directed graph
  in: in degree == 1
  out: out degree == 2
  height: vertex height (time units)
  seq: nucleotide sequence (if known), NaN otherwise
  """
  in::Bool
  out::Bool
  height::Float64
  seq
end

TreeVertex(h::Float64, s::Nucleotide2bitSeq) = TreeVertex(true, false, h, s)

function display(object::TreeVertex)
  """
  Display for TreeVertex
  """
  if object.in == true
    if object.out == false
      println("Leaf @ $(object.height)")
    else
      println("Node @ $(object.height)")
    end
  else
    if object.out == true
      println("Root @ $(object.height)")
    else
      error("Invalid TreeVertex")
    end
  end
end

function display(object::Vector{TreeVertex})
  """
  Display for TreeVertex
  """
  for i = 1:length(object)
    if object[i].in == true
      if object[i].out == false
        print("Leaf ")
      else
        print("Node ")
      end
    else
      if object[i].out == true
        print("Root ")
      else
        error("Invalid TreeVertex")
      end
    end
  end
end

function push!(a::TreeVertex, b::TreeVertex)
  a = [a, b]
end

type TreeEdge
  """
  Describes an edge in a phylogenetic tree directed graph
  from: vertex identification
  to: vertex identification
  """
  from::Int64
  to::Int64
end

type Tree
  """
  All possible phylogenetic trees described by a directed graph
  """
  Vertices::Vector{TreeVertex}
  Edges::Vector{TreeEdge}
end

type SEIR_priors{T<:UnivariateDistribution}
  """
  Prior distributions for SEIR model inference
  """
  α::T
  β::T
  ρ::T
  γ::T
  η::T
  ν::T
end

type SEIR_actuality
  """
  Contains actual event times
  """
  exposed::Vector{Float64}
  infectious::Vector{Float64}
  removed::Vector{Float64}
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
  infectious::Vector{Float64}
  exposed::Vector{Float64}
  removed::Vector{Float64}
end

type SEIR_trace
  """
  Contains an MCMC trace object
  """
  α::Vector{Float64}
  β::Vector{Float64}
  ρ::Vector{Float64}
  γ::Vector{Float64}
  η::Vector{Float64}
  ν::Vector{Float64}
  aug::Vector{SEIR_augmented}
  sources::Vector{Array{Float64}}
  logposterior::Vector{Float64}
end
