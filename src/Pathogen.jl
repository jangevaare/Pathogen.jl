module Pathogen

using Distributions, DataFrames

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

include("Types.jl")
include("Sequence.jl")
include("Substitution_models.jl")

end
