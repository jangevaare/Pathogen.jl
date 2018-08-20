function _plot_pathway(pop::Population,
                       network::TransmissionNetwork,
                       i::Int64)
  if !network.external[i] && any(network.internal[:, i])
    source = findfirst(network.internal[:, i])
    x = pop.risks[:x][[i; source]]
    y = pop.risks[:y][[i; source]]
    return x, y
  else
    return Float64[], Float64[]
  end
end
