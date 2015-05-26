"""
utilities.jl - basic utilities for dealing with sequence data
Justin Angevaare
May 2015
"""

function generate_sequence(n::Int, π_A::Float64, π_G::Float64, π_C::Float64, π_U::Float64)
"""
Generate a nucleotide sequence of length `n`, with specific nucleotide frequencies
"""
  @assert(sum([π_A, π_G, π_C, π_U]) == 1, "Nucleotide frequencies must sum to 1")
  @assert(all(0 .< [π_A, π_G, π_C, π_U] .< 1), "Each nucleotide frequency must be between 0 and 10")
  sequence = fill(0, n)
  for i = 1:n
    sequence[i]=find(rand(Multinomial(1, [π_A, π_G, π_C, π_U])))[1]
  end
  return sequence
end
