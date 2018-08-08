function _pathway_from(id::Int64,
                       network::TransmissionNetwork;
                       depth::Real=Inf,
                       debug_level::Int64=0)
  path = Int64[]
  if network.external[id] | any(network.internal[:, id])
    push!(path, id)
    path_lengths = [0]
    while (length(path) > path_lengths[end]) & (depth >= length(path_lengths))
      push!(path_lengths, length(path))
      for j in path[(path_lengths[end-1]+1):path_lengths[end]]
        append!(path, find(network.internal[j, :]))
      end
    end
  end
  if debug_level >= 3
    println("_pathway_from: transmission pathway from $id: $path")
  end
  return path
end
