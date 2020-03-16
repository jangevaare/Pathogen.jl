function TransmissionNetworkDistribution(x::MarkovChain)
  return TNDistribution(x.transmission_network)
end

function TransmissionNetworkDistribution(iter, x::MarkovChain)
  return TNDistribution(x.transmission_network[iter])
end

function TransmissionNetworkDistribution(x::MCMC)
  return TNDistribution([TNDistribution(y) for y in x.markov_chains])
end
