function _plot_pathway(pop::Population,
                       network::TransmissionNetwork,
                       i::Int64)
  if !network.external[i] && any(network.internal[:, i])
    source = findfirst(network.internal[:, i])
    x = pop.risks[[i; source], :x]
    y = pop.risks[[i; source], :y]
    return x, y
  else
    return Float64[], Float64[]
  end
end
