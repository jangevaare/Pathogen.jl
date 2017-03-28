"""
pathwayto(individual::Int64,
          network::Network,
          depth=Inf::Real)

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
