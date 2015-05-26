module Pathogen

using Distributions, DataFrames, BioSeq

export
  # Types.jl
  population,

  # Sequence.jl
  nucleotide_convert,
  nucleotide_revert,
  generate_sequence,
  mutate_sequence,

  # Substitution_models.jl
  JC69,
  K80

include("types.jl")
include("utilites.jl")
include("substitute.jl")
include("simulate.jl")

end
