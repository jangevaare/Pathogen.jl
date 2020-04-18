function _pathway_from(id::I,
                       network::TransmissionNetwork;
                       depth::Real=Inf) where {I <: Integer}
  path = I[]
  if network.external[id] | any(network.internal[:, id])
    push!(path, id)
    path_lengths = [0]
    while (length(path) > path_lengths[end]) & (depth >= length(path_lengths))
      push!(path_lengths, length(path))
      for j in path[(path_lengths[end-1]+1):path_lengths[end]]
        append!(path, findall(network.internal[j, :]))
      end
    end
  end
  @debug "_pathway_from: transmission pathway from $id: $path"
  return path
end
