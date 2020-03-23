struct TransmissionNetworkDistribution
  external::Array{Float64, 1}
  internal::Array{Float64, 2}

  function TransmissionNetworkDistribution(x::Vector{TransmissionNetwork})
    return new(mean([y.external for y in x]), mean([y.internal for y in x]))
  end

  function TransmissionNetworkDistribution(x::Vector{TransmissionNetworkDistribution})
    return new(mean([y.external for y in x]), mean([y.internal for y in x]))
  end
end

const TNDistribution = TransmissionNetworkDistribution
const TransmissionNetworkPrior = TransmissionNetworkDistribution
const TNPrior = TransmissionNetworkPrior
const TransmissionNetworkPosterior = TransmissionNetworkDistribution
const TNPosterior = TransmissionNetworkPosterior

function Base.sum(x::TNDistribution)
  return sum(x.external) + sum(x.internal)
end

function Base.show(io::IO, x::TNDistribution)
  return print(io, "Transmission network distribution")
end
