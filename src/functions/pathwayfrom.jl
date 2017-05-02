"""
pathwayfrom(individual::Int64,
            network::Network,
            depth=Inf::Real)

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
