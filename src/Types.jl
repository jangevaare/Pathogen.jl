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

type SEIR_events
  """
  Contains observed and actual event times
  """
  ids::Vector{Int64}
  exposed_actual::Vector{Float64}
  infectious_actual::Vector{Float64}
  infectious_observed::Vector{Float64}
  removed_actual::Vector{Float64}
  removed_observed::Vector{Float64}
  covariates = Vector{Vector{Any}}
  seq = Vector{Any}
end

type SEIR_augmented
  """
  Contains event times from data augmentation
  """
  infectious_augmented::Vector{Float64}
  exposed_augmented::Vector{Float64}
  recovered_augmented::Vector{Float64}
end
