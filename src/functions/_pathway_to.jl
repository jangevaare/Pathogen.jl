function _pathway_to(id::Int64,
                     network::TransmissionNetwork;
                     depth::Real=Inf,
                     debug_level::Int64=0)
  path = Int64[]
  if network.external[id] | any(network.internal[:, id])
    push!(path, id)
    next_id = findfirst(network.internal[:, id])
    while (length(path) <= depth) & (next_id != 0)
      push!(path, next_id)
      next_id = findfirst(network.internal[:, path[end]])
    end
  end
  if debug_level >= 3
    println("_pathway_to: transmission pathway to $id: $path")
  end
  return path
end
