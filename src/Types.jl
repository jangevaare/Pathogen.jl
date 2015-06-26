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

type Tree
  """
  Phylogenetic tree
  """
  sequences::Vector{Nucleotide2bitSeq}
  positions::Vector{Vector{Bool}}
  distances::Vector{Vector{Float64}}
end
