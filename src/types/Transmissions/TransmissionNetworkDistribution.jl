struct TransmissionNetworkDistribution
  external::Array{Float64, 1}
  internal::Array{Float64, 2}

  function TransmissionNetworkDistribution(x::Vector{TransmissionNetwork})
    return new(mean([y.external for y in x]), mean([y.internal for y in x]))
  end
end

const TNDistribution = TransmissionNetworkDistribution
const TransmissionNetworkPrior = TransmissionNetworkDistribution
const TNPrior = TransmissionNetworkPrior
const TransmissionNetworkPosterior = TransmissionNetworkDistribution
const TNPosterior = TransmissionNetworkPosterior

function TransmissionNetworkDistribution(x::MarkovChain)
  return TNDistribution(x.transmission_network)
end

function TransmissionNetworkDistribution(iter, x::MarkovChain)
  return TNDistribution(x.transmission_network[iter])
end

function TransmissionNetworkDistribution(x::MCMC)
  return TNDistribution([TNDistribution(y) for y in x.markov_chains])
end

function Base.sum(x::TNDistribution)
  return sum(x.external) + sum(x.internal)
end

function Base.show(io::IO, x::TNDistribution)
  return print(io, "Transmission network distribution")
end
