"""
Return the transmission pathway leading to an individual
"""
function pathwayto(individual::Int64,
                   network::Network,
                   depth=Inf::Real)
  path = Int64[]
  nextindividual = findfirst(network.internal[:, individual])
  if network.external[individual]
    push!(path, individual)
  elseif nextindividual != 0
    push!(path, individual)
    while length(path) <= depth && nextindividual != 0
      push!(path, nextindividual)
      nextindividual = findfirst(network.internal[:, path[end]])
    end
  end
  return path
end


"""
Return the transmission pathway leading from an individual
"""
function pathwayfrom(individual::Int64,
                     network::Network,
                     depth=Inf::Real)
   path = Int64[]
   if network.external[individual] | any(network.internal[:, individual])
     push!(path, individual)
     pathlengths = [0]
     while (length(path) > pathlengths[end]) & (depth >= length(pathlengths))
       push!(pathlengths, length(path))
       for j in path[(pathlengths[end-1]+1):pathlengths[end]]
         append!(path, find(network.internal[j, :]))
       end
     end
   end
  return path
end
