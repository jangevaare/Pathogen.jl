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
  height::FloatingPoint
  seq
end

TreeVertex(i::Bool, o::Bool) = TreeVertex(i, o, 0., NaN)

TreeVertex(i::Bool, o::Bool, h::FloatingPoint) = TreeVertex(i, o, h, NaN)

TreeVertex(i::Bool, o::Bool, h::FloatingPoint) = TreeVertex(i, o, h, NaN)

TreeVertex() = TreeVertex(false, true, 0., NaN)

TreeVertex(h::FloatingPoint) = TreeVertex(true, true, h, NaN)

TreeVertex(h::FloatingPoint, s::Nucleotide2bitSeq) = TreeVertex(true, false, h, s)

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

