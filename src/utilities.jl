import Base.convert


import Base.maximum


"""
Add a conversion method to move from Nucleotide to an integer vector
"""
function convert(::Type{Vector{Int64}}, x::Nucleotide2bitSeq)
  return sub2ind((2,2), x.b1 .+1, x.b2 .+1)
end


"""
Add a conversion method to move from an integer vector to a nucleotide sequence
"""
function convert(::Type{Nucleotide2bitSeq}, x::Vector{Int64})
  b1,b2 = ind2sub((2,2), x)
  if length(x) == 1
    return Nucleotide2bitSeq(convert(BitArray, [b1 - 1]), convert(BitArray, [b2 - 1]))
  else
    return Nucleotide2bitSeq(convert(BitArray, b1 - 1), convert(BitArray, b2 - 1))
  end
end


"""
Add a conversion method to move from nucleotide base to an integer
"""
function convert(::Type{Int64}, x::Nucleotide2bitBase)
  return sub2ind((2,2), x.b1 .+1, x.b2 .+1)
end


"""
Add a conversion method to move from an integer to a nucleotide base
"""
function convert(::Type{Nucleotide2bitBase}, x::Int64)
  b1,b2 = ind2sub((2,2), x)
  return Nucleotide2bitBase(convert(Bool, b1 - 1), convert(Bool, b2 - 1))
end


"""
Find the disease state of a specific individual
"""
function findstate(population::Population, individual::Int64, time::Float64)
  if sum(population.events[individual][1] .< time) == 0
    return "S"
  elseif sum(population.events[individual][1] .< time) > sum(population.events[individual][3] .< time)
    return "E"
  elseif sum(population.events[individual][3] .< time) > sum(population.events[individual][4] .< time)
    return "I"
  elseif sum(population.events[individual][1] .< time) == sum(population.events[individual][4] .< time)
    return "R"
  else
    return "unknown"
  end
end


"""
Find the disease state of a specific individual
"""
function findstate(trace::SEIR_trace, iteration::Int64, individual::Int64, time::Float64)
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


"""
Find the disease state of a specific individual
"""
function findstate(trace::SIR_trace, iteration::Int64, individual::Int64, time::Float64)
  if isnan(trace.aug[iteration].infectious[individual]) || trace.aug[iteration].infectious[individual] > time
    return "S"
  elseif trace.aug[iteration].infectious[individual] < time && isnan(trace.aug[iteration].removed[individual]) || trace.aug[iteration].removed[individual] > time
    return "I"
  elseif trace.aug[iteration].removed[individual] < time
    return "R"
  else
    return "unknown"
  end
end


"""
Find the disease state of a specific individual
"""
function findstate(actual::SEIR_actual, individual::Int64, time::Float64)
  if isnan(actual.exposed[individual]) || actual.exposed[individual] > time
    return "S"
  elseif actual.exposed[individual] < time && isnan(actual.infectious[individual]) || actual.infectious[individual] > time
    return "E"
  elseif actual.infectious[individual] < time && isnan(actual.removed[individual]) || actual.removed[individual] > time
    return "I"
  elseif actual.removed[individual] < time
    return "R"
  else
    return "unknown"
  end
end


"""
Find the disease state of a specific individual
"""
function findstate(actual::SIR_actual, individual::Int64, time::Float64)
  if isnan(actual.infectious[individual]) || actual.infectious[individual] > time
    return "S"
  elseif actual.infectious[individual] < time && isnan(actual.removed[individual]) || actual.removed[individual] > time
    return "I"
  elseif actual.removed[individual] < time
    return "R"
  else
    return "unknown"
  end
end


"""
Create dataframes with all necessary plotting information
"""
function plotdata(population::Population, time::Float64)
  states = DataFrame(id = fill(NaN,4),
                     x = fill(NaN,4),
                     y = fill(NaN,4),
                     state = ["S", "E", "I", "R"])
  routes = DataFrame(x = Float64[],
                     y = Float64[],
                     age = Float64[])
  for i = 2:length(population.events)
    states = vcat(states, DataFrame(x = population.history[i][1][find(time .> population.events[i][5])[end]][1],
                                    y = population.history[i][1][find(time .> population.events[i][5])[end]][2],
                                    state = findstate(population, i, time)))
    for j = 1:sum(population.events[i][1].< time)
      source = population.events[i][2][j]
      age = time - population.events[i][1][j]
      if source > 1
        routes = vcat(routes, DataFrame(x = population.history[i][1][find(time .> population.events[i][5])[end]][1],
                                        y = population.history[i][1][find(time .> population.events[i][5])[end]][2],
                                        age = "$age"))
        routes = vcat(routes, DataFrame(x = population.history[source][1][find(time .> population.events[source][5])[end]][1],
                                        y = population.history[source][1][find(time .> population.events[source][5])[end]][2],
                                        age = "$age"))
      end
    end
  end
  return states, routes
end


"""
Create dataframes with all necessary plotting information
"""
function plotdata(obs::SEIR_observed,
                  trace::SEIR_trace,
                  iteration::Int64,
                  time::Float64)
  states = DataFrame(id = fill(NaN,4),
                     x = fill(NaN,4),
                     y = fill(NaN,4),
                     state = ["S", "E", "I", "R"])
  routes = DataFrame(x = Float64[],
                     y = Float64[],
                     age = Float64[])
  for i = 1:length(obs.infectious)
    states = vcat(states, DataFrame(x = obs.covariates[i][1],
                                    y = obs.covariates[i][2],
                                    state = findstate(trace, iteration, i, time)))
    if states[:state][end] != "S"
      source = findfirst(trace.network[iteration][:,i])-1
      age = time - trace.aug[iteration].exposed[i]
        if source > 1
        routes = vcat(routes, DataFrame(x = obs.covariates[i][1],
                                        y = obs.covariates[i][2],
                                        age = "$age"))
        routes = vcat(routes, DataFrame(x = obs.covariates[source][1],
                                        y = obs.covariates[source][2],
                                        age = "$age"))
      end
    end
  end
  return states, routes
end


"""
Create dataframes with all necessary plotting information
"""
function plotdata(obs::SIR_observed,
                  trace::SIR_trace,
                  iteration::Int64,
                  time::Float64)
  states = DataFrame(id = fill(NaN,3),
                     x = fill(NaN,3),
                     y = fill(NaN,3),
                     state = ["S", "I", "R"])
  routes = DataFrame(x = Float64[],
                     y = Float64[],
                     age = Float64[])
  for i = 1:length(obs.infectious)
    states = vcat(states, DataFrame(x = obs.covariates[i][1],
                                    y = obs.covariates[i][2],
                                    state = findstate(trace, iteration, i, time)))
    if states[:state][end] != "S"
      source = findfirst(trace.network[iteration][:,i])-1
      age = time - trace.aug[iteration].infectious[i]
        if source > 1
        routes = vcat(routes, DataFrame(x = obs.covariates[i][1],
                                        y = obs.covariates[i][2],
                                        age = "$age"))
        routes = vcat(routes, DataFrame(x = obs.covariates[source][1],
                                        y = obs.covariates[source][2],
                                        age = "$age"))
      end
    end
  end
  return states, routes
end


"""
Create dataframes with all necessary plotting information
"""
function plotdata(actual::SEIR_actual,
                  population::Population,
                  time::Float64)
  states = DataFrame(id = fill(NaN, 4),
                     x = fill(NaN, 4),
                     y = fill(NaN, 4),
                     state = ["S", "E", "I", "R"])
  routes = DataFrame(x = Float64[],
                     y = Float64[],
                     age = Float64[])
  for i = 1:length(actual.infectious)
    states = vcat(states, DataFrame(x = actual.covariates[i][1],
                                    y = actual.covariates[i][2],
                                    state = findstate(actual, i, time)))
    if states[:state][end] != "S"
      source = population.events[i][2][j]
      age = time - actual.exposed[i]
        if source > 1
        routes = vcat(routes, DataFrame(x = actual.covariates[i][1],
                                        y = actual.covariates[i][2],
                                        age = "$age"))
        routes = vcat(routes, DataFrame(x = actual.covariates[source][1],
                                        y = actual.covariates[source][2],
                                        age = "$age"))
      end
    end
  end
  return states, routes
end


"""
Create dataframes with all necessary plotting information
"""
function plotdata(actual::SIR_actual,
                  population::Population,
                  time::Float64)
  states = DataFrame(id = fill(NaN, 3),
                     x = fill(NaN, 3),
                     y = fill(NaN, 3),
                     state = ["S", "I", "R"])
  routes = DataFrame(x = Float64[],
                     y = Float64[],
                     age = Float64[])
  for i = 1:length(actual.infectious)
    states = vcat(states, DataFrame(x = actual.covariates[i][1],
                                    y = actual.covariates[i][2],
                                    state = findstate(actual, i, time)))
    if states[:state][end] != "S"
      source = population.events[i][2][j]
      age = time - actual.exposed[i]
        if source > 1
        routes = vcat(routes, DataFrame(x = actual.covariates[i][1],
                                        y = actual.covariates[i][2],
                                        age = "$age"))
        routes = vcat(routes, DataFrame(x = actual.covariates[source][1],
                                        y = actual.covariates[source][2],
                                        age = "$age"))
      end
    end
  end
  return states, routes
end


"""
Find maximum augmented event time
"""
function maximum(aug::SEIR_augmented)
  return maximum([aug.exposed aug.infectious aug.removed])
end


"""
Find maximum augmented event time
"""
function maximum(aug::SIR_augmented)
  return maximum([aug.infectious aug.removed])
end


"""
Find maximum observed event time
"""
function maximum(obs::Observed)
  return maximum([obs.infectious obs.removed])
end


"""
Return the transmission pathways leading to an individual
"""
function pathwayto(infection::Int64, network::Array{Bool,2}, debug=false::Bool)
  path = [infection]
  while path[end] > 0
    push!(path, findfirst(network[:, path[end]])-1)
  end
  if debug
    println("Pathway to: $path")
  end
  return path
end


"""
Return all transmission pathways leading to specified individuals
"""
function pathwaysto(infections::Vector{Int64}, network::Array{Bool,2}, debug=false::Bool)
  paths = Array[Int64[]]
  for i = 1:length(infections)
    if i > 1
      push!(paths, Int64[])
    end
    push!(paths[i], infections[i])
    while paths[i][end] > 0
      push!(paths[i], findfirst(network[:, paths[i][end]])-1)
    end
    if debug
      println("Pathway to: $(paths[i])")
    end
  end
  return paths
end


"""
Return all transmission pathways
"""
function pathwaysto(network::Array{Bool,2}, debug=false::Bool)
  infections = find(sum(network,1))
  return pathwaysto(infections, network, debug)
end


"""
Return the transmission pathways leading from an individual
"""
function pathwayfrom(infection::Int64, network::Array{Bool,2}, depth=0::Int64, debug=false::Bool)
  path = [infection]
  pathlengths = [0]
  if depth == 0
    while length(path) > pathlengths[end]
      push!(pathlengths, length(path))
      for j in path[(pathlengths[end-1]+1):pathlengths[end]]
        append!(path, find(network[j+1,:]))
      end
    end
  else
    while depth + 1 > length(path) > pathlengths[end]
      push!(pathlengths, length(path))
      for j in path[(pathlengths[end-1]+1):pathlengths[end]]
        append!(path, find(network[j+1,:]))
      end
    end
  end

  if debug
    println("Pathway from: $path")
  end
  return path
end


pathwayfrom(infection::Int64, network::Array{Bool,2}, debug=false::Bool) = pathwayfrom(infection, network, 0, debug)


"""
Return all transmission pathways leading from specified individuals
"""
function pathwaysfrom(infections::Vector{Int64}, network::Array{Bool,2}, depth=0::Int64, debug=false::Bool)
  paths = Array[Int64[]]
  for i = 1:length(infections)
    if i > 1
      push!(paths, Int64[])
    end
    push!(paths[i], infections[i])
    pathlengths = [0]
    if depth == 0
      while length(paths[i]) > pathlengths[end]
        push!(pathlengths, length(paths[i]))
        for j in paths[i][(pathlengths[end-1]+1):pathlengths[end]]
          append!(paths[i], find(network[j+1,:]))
        end
      end
    else
      while depth + 1 > length(paths[i]) > pathlengths[end]
        push!(pathlengths, length(paths[i]))
        for j in paths[i][(pathlengths[end-1]+1):pathlengths[end]]
          append!(paths[i], find(network[j+1,:]))
        end
      end
    end
    if debug
      println("Pathway from: $(paths[i])")
    end
  end
  return paths
end


"""
Return all transmission pathways
"""
function pathwaysfrom(network::Array{Bool,2}, depth=0::Int64, debug=false::Bool)
  infections = find(sum(network,1))
  return pathwaysfrom(infections, network, depth, debug)
end

pathwaysfrom(network::Array{Bool,2}, debug::Bool) = pathwaysfrom(network, 0, debug)
