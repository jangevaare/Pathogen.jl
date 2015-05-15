"""
Sequence.jl - basic utilities for dealing with sequence data
Justin Angevaare
May 2015
"""

function nucleotide_convert(x::ASCIIString)
"""
Convert a sequence of nucleotide bases into an integer representation
"""
  sequence = fill(0, length(x))
  for i = 1:length(x)
    if x[i] == 'A' sequence[i] = 1 end
    if x[i] == 'G' sequence[i] = 2 end
    if x[i] == 'C' sequence[i] = 3 end
    if x[i] == 'U' sequence[i] = 4 end
  end
  return sequence
end

function nucleotide_revert(x::Vector{Int})
"""
Convert an integer sequence into a nucleotide base representation
"""
  sequence = fill("Z", length(x))
  for i = 1:length(x)
    if x[i] == 1 sequence[i] = "A" end
    if x[i] == 2 sequence[i] = "G" end
    if x[i] == 3 sequence[i] = "C" end
    if x[i] == 4 sequence[i] = "U" end
  end
sequence
end

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
