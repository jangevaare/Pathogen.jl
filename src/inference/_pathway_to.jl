function _pathway_to(id::I,
                     network::TransmissionNetwork;
                     depth::Real=Inf) where {I <: Integer}
  path = I[]
  if network.external[id] | any(network.internal[:, id])
    push!(path, id)
    next_id = findfirst(network.internal[:, id])
    while (length(path) <= depth) & (typeof(next_id) != Nothing)
      push!(path, next_id)
      next_id = findfirst(network.internal[:, path[end]])
    end
  end
  @debug "_pathway_to: transmission pathway to $id: $path"
  return path
end
