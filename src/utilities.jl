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
  @assert(all(0 .< [π_A, π_G, π_C, π_U] .< 1), "Each nucleotide frequency must be between 0 and 1")
  sequence = fill(0, n)
  for i = 1:n
    sequence[i]=findfirst(rand(Multinomial(1, [π_A, π_G, π_C, π_U])))
  end
  return nucleotide("AGCU"[sequence])
end

function find_state(population::Population, individual::Int64, time::Float64)
  """
  Find the disease state of a specific individual
  """
  # exposure times, exposure source, infection times, recovery times, covariate times, sequence times
  if length(population.events[individual][1]) == 0
    return "S"
  elseif length(population.events[individual][1]) > length(population.events[individual][2])
    return "E"
  elseif length(population.events[individual][2]) > length(population.events[individual][3])
    return "I"
  elseif length(population.events[individual][1]) == length(population.events[individual][4])
    return "S*"
  else
    return "unknown"
  end
end
