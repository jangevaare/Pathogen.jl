"""
Return the transmission pathway leading to an individual
"""
function pathwayto(individual::Int64,
                   network::Network,
                   depth=Inf::Real)
  path = Int64[]
  nextindividual = findfirst(network.external[:, individual])
  if network.internal[individual]
    push!(path, individual)
  elseif nextindividual != 0
    push!(path, individual)
    while length(path) <= depth && nextindividual != 0
      push!(path, nextindividual)
      nextindividual = findfirst(network.external[:, path[end]])
    end
  end
  return path
end


"""
Return all transmission pathways leading to specified individuals
"""
function pathwayto(individuals::Vector{Int64},
                   network::Network,
                   depth=Inf::Real)
  paths = Array[Int64[]]
  for i in individuals
    push!(paths, pathwayto(individual, network, depth))
  end
  return paths
end


"""
Return all transmission pathways leading to individuals
"""
function pathwayto(network::Network,
                   depth=Inf::Real)
  return pathwayto(1:length(network.internal), network, depth)
end


"""
Return the transmission pathway leading from an individual
"""
function pathwayfrom(individual::Int64,
                     network::Network,
                     depth=Inf::Real)
   path = Int64[]
   push!(path, individual)
   pathlengths = [0]
   while depth >= length(path) > pathlengths[end]
     push!(pathlengths, length(path))
     for j in path[(pathlengths[end-1]+1):pathlengths[end]]
       append!(path, find(network.external[j, :]))
     end
   end
  return path
end


"""
Return all transmission pathways leading from specified individuals
"""
function pathwayfrom(individuals::Vector{Int64},
                     network::Network,
                     depth=Inf::Real)
  paths = Array[Int64[]]
  for i in individuals
    push!(paths, pathwayfrom(i, network, depth))
  end
  return paths
end


"""
Return all transmission pathways leading from individuals
"""
function pathwayfrom(network::Network,
                     depth=Inf::Real)
  return pathwayfrom(1:length(network.internal), network, depth)
end
