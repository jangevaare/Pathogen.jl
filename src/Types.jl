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

TreeVertex() = TreeVertex(false, true, 0., NaN)

TreeVertex(h::Float64) = TreeVertex(true, true, h, NaN)

TreeVertex(h::Float64, s::Nucleotide2bitSeq) = TreeVertex(true, false, h, s)

TreeVertex(i::Bool, o::Bool) = TreeVertex(i, o, 0., NaN)

TreeVertex(i::Bool, o::Bool, h::Float64) = TreeVertex(i, o, h, NaN)

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

