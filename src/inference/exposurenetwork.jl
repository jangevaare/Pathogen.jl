"""
Probablistically generate a network object based on exposure network rates
"""
function rand(network_rates::Array{Array{Float64}, 1})
  external_rates = network_rates[1]
  internal_rates = network_rates[2]
  external_network = fill(false, length(external_rates))
  internal_network = fill(false, size(internal_rates))
  if !(length(external_rates) == size(internal_rates, 1) == size(internal_rates, 2))
    throw(BoundsError)
  end
  for i = 1:length(external_rates)
    external_total = external_rates[i]
    internal_total = sum(internal_rates[:, i])
    if rand() < external_total/(external_total + internal_total)
      external_network[i] = true
    else
      source = rand(Multinomial(1, network_rates.internal[:, i]/internal_total))
      internal_network[source, i] = true
    end
  end
  return Network(external_network, internal_network)
end


"""
Propose an exposure network based on a previous exposure network and exposure
network rates
"""
function propose(individuals::Vector{Int64},
                 network::Network,
                 network_rates::Array{Array{Float64}, 1})
  proposal = network
  for i in individuals
    proposal.external[i] = false
    proposal.internal[:, i] = false
    external_total = network_rates.external[i]
    internal_total = sum(network_rates.internal[:, i])
    if rand() < external_total/(external_total + internal_total)
      proposal.external[i] = true
    else
      source = rand(Multinomial(1, network_rates.internal[:, i]/internal_total))
      proposal.internal[source, i] = true
    end
  end
  return proposal
end
