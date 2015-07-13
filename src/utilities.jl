"""
utilities.jl - basic utilities required by Pathogen.jl
Justin Angevaare
May 2015
"""
import Base.convert

function convert(::Type{Vector{Int64}}, x::Nucleotide2bitSeq)
  """
  Add a conversion method to move from Nucleotide to an integer vector
  """
  return sub2ind((2,2), x.b1 .+1, x.b2 .+1)
end

function convert(::Type{Nucleotide2bitSeq}, x::Vector{Int64})
  """
  Add a conversion method to move from an integer vector to a nucleotide sequence
  """
  b1,b2 = ind2sub((2,2), x)
  if length(x) == 1
    return Nucleotide2bitSeq(convert(BitArray, [b1 - 1]), convert(BitArray, [b2 - 1]))
  else
    return Nucleotide2bitSeq(convert(BitArray, b1 - 1), convert(BitArray, b2 - 1))
  end
end

function convert(::Type{Int64}, x::Nucleotide2bitBase)
  """
  Add a conversion method to move from nucleotide base to an integer
  """
  return sub2ind((2,2), x.b1 .+1, x.b2 .+1)
end

function convert(::Type{Nucleotide2bitBase}, x::Int64)
  """
  Add a conversion method to move from an integer to a nucleotide base
  """
  b1,b2 = ind2sub((2,2), x)
  return Nucleotide2bitBase(convert(Bool, b1 - 1), convert(Bool, b2 - 1))
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

function plotdata(population::Population, time::Float64)
  """
  Create dataframes with all necessary plotting information
  """
  states = DataFrame(id = fill(NaN,4), x = fill(NaN,4), y = fill(NaN,4), state = ["S", "E", "I", "S*"])
  routes = DataFrame(x = Float64[], y = Float64[], age = Float64[])
  for i = 2:length(population.events)
    states = vcat(states, DataFrame(x = population.history[i][1][1,find(time .> population.events[i][5])[end]], y = population.history[i][1][2,find(time .> population.events[i][5])[end]], state = findstate(population, i, time)))
    for j = 1:sum(population.events[i][1].< time)
      source = population.events[i][2][j]
      age = time - population.events[i][1][j]
      if source > 1
        routes = vcat(routes, DataFrame(x = population.history[i][1][1,find(time .> population.events[i][5])[end]], y = population.history[i][1][2,find(time .> population.events[i][5])[end]], age = "$age"))
        routes = vcat(routes, DataFrame(x = population.history[source][1][1,find(time .> population.events[source][5])[end]], y = population.history[source][1][2,find(time .> population.events[source][5])[end]], age = "$age"))
      end
    end
  end
  return states, routes
end
