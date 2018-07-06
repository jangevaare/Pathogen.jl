function _pathway_to(population::DataFrame,
                     network::TransmissionNetwork,
                     i::Int64)
  if !network.external[i] && any(network.internal[:, i])
    source = findfirst(network.internal[:, i])
    x = population[:x][[i; source]]
    y = population[:y][[i; source]]
    return x, y
  else
    return Float64[], Float64[]
  end
end
