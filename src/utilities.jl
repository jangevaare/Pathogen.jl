"""
utilities.jl
"""

import Base.convert
import Base.maximum

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

function findstate(trace::SEIR_trace, iteration::Int64, individual::Int64, time::Float64)
  """
  Find the disease state of a specific individual
  """
  if isnan(trace.aug[iteration].exposed[individual]) || trace.aug[iteration].exposed[individual] > time
    return "S"
  elseif trace.aug[iteration].exposed[individual] < time && isnan(trace.aug[iteration].infectious[individual]) || trace.aug[iteration].infectious[individual] > time
    return "E"
  elseif trace.aug[iteration].infectious[individual] < time && isnan(trace.aug[iteration].removed[individual]) || trace.aug[iteration].removed[individual] > time
    return "I"
  elseif trace.aug[iteration].removed[individual] < time
    return "R"
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
    states = vcat(states, DataFrame(x = population.history[i][1][find(time .> population.events[i][5])[end]][1], y = population.history[i][1][find(time .> population.events[i][5])[end]][2], state = findstate(population, i, time)))
    for j = 1:sum(population.events[i][1].< time)
      source = population.events[i][2][j]
      age = time - population.events[i][1][j]
      if source > 1
        routes = vcat(routes, DataFrame(x = population.history[i][1][find(time .> population.events[i][5])[end]][1], y = population.history[i][1][find(time .> population.events[i][5])[end]][2], age = "$age"))
        routes = vcat(routes, DataFrame(x = population.history[source][1][find(time .> population.events[source][5])[end]][1], y = population.history[source][1][find(time .> population.events[source][5])[end]][2], age = "$age"))
      end
    end
  end
  return states, routes
end

function plotdata(obs::SEIR_observed, trace::SEIR_trace, iteration::Int64, time::Float64)
  """
  Create dataframes with all necessary plotting information
  """
  states = DataFrame(id = fill(NaN,4), x = fill(NaN,4), y = fill(NaN,4), state = ["S", "E", "I", "R"])
  routes = DataFrame(x = Float64[], y = Float64[], age = Float64[])
  for i = 1:length(obs.infectious)
    states = vcat(states, DataFrame(x = obs.covariates[i][1], y = obs.covariates[i][2], state = findstate(trace::SEIR_trace, iteration, i, time)))
    if states[:state][end] != "S"
      source = findfirst(trace.network[iteration][:,i])-1
      age = time - trace.aug[iteration].exposed[i]
        if source > 1
        routes = vcat(routes, DataFrame(x = obs.covariates[i][1], y = obs.covariates[i][2], age = "$age"))
        routes = vcat(routes, DataFrame(x = obs.covariates[source][1], y = obs.covariates[source][2], age = "$age"))
      end
    end
  end
  return states, routes
end

function maximum(aug::SEIR_augmented)
  """
  Find maximum augmented event time
  """
  return maximum([aug.exposed aug.infectious aug.removed])
end

# function isseq(x::Vector{Any})
#   """
#   Identifies the pressence of a Nucleotide2bitSeq in a Vector{Any}
#   """
#   seq = fill(false, length(x))
#   for i = 1:length(x)
#     if typeof(x[i]) == Nucleotide2bitSeq
#       seq[i] = true
#     end
#   end
#   return seq
# end

function isseq(x::Vector{Any})
  """
  Identifies the pressence of a Vector{Int64} in a Vector{Any}
  """
  seq = fill(false, length(x))
  for i = 1:length(x)
    if typeof(x[i]) == Vector{Int64}
      seq[i] = true
    end
  end
  return seq
end
