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

function findstate(population::Population, individual::Int64, time::Float64)
  """
  Find the disease state of a specific individual
  """
  if sum(population.events[individual][1] .< time) == 0
    return "S"
  elseif sum(population.events[individual][1] .< time) > sum(population.events[individual][3] .< time)
    return "E"
  elseif sum(population.events[individual][3] .< time) > sum(population.events[individual][4] .< time)
    return "I"
  elseif sum(population.events[individual][1] .< time) == sum(population.events[individual][4] .< time)
    return "S*"
  else
    return "unknown"
  end
end

function plotdata(population, time)
  """
  Create dataframes with all necessary plotting information
  """
#   states = DataFrame(x = Float64[], y = Float64[], state = ASCIIString[])
  states = DataFrame(id = fill(NaN,4), x = fill(NaN,4), y = fill(NaN,4), state = ["S", "E", "I", "S*"])
  routes = DataFrame(x = Float64[], y = Float64[], age = Float64[])
  for i = 2:length(population.events)
    states = vcat(states, DataFrame(x = population.history[i][1][1,(time .> population.events[i][5])[end]], y = population.history[i][1][2,(time .> population.events[i][5])[end]], state = findstate(population, i, time)))
    for j = 1:sum(population.events[i][1].< time)
      source = population.events[i][2][j]
      age = time - population.events[i][1][j]
      if source > 1
        routes = vcat(routes, DataFrame(x = population.history[i][1][1,(time .> population.events[i][5])[end]], y = population.history[i][1][2,(time .> population.events[i][5])[end]], age = "$age"))
        routes = vcat(routes, DataFrame(x = population.history[source][1][1,(time .> population.events[source][5])[end]], y = population.history[source][1][2,(time .> population.events[source][5])[end]], age = "$age"))
      end
    end
  end
  return states, routes
end
