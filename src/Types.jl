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

# abstract TreeFeature

# type TreeNode <: TreeFeature
#   """
#   Node in a phylogenetic tree
#   `distance` is from ancestral node
#   `branches` are further nodes or leaves (vector of length 2 for bifurcation)
#   """
#   distance::Float64
#   branches::Vector{TreeFeature}
# end

# type TreeLeaf <: TreeFeature
#   """
#   Leaf in a phylogenetic tree
#   `distance` is from ancenstral node
#   `id` is a single sequence identifier
#   """
#   distance::Float64
#   id::Int64
# end

# type Tree
#   """
#   Phylogenetic tree
#   """
#   structure::TreeNode
#   sequences::Vector{Nucleotide2bitSeq}
# end

type Tree
  """
  Phylogenetic tree
  """
  sequence::Vector{Nucleotide2bitSeq}
  location::Vector{Vector{Bool}}
  distance::Vector{Vector{Float64}}
end
