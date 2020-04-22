struct TransmissionNetworkDistribution <: AbstractTransmissionNetwork
  external::Array{Float64, 1}
  internal::Array{Float64, 2}

  function TransmissionNetworkDistribution(e::Vector{Float64},
                                           i::Array{Float64, 2})
    @boundscheck if !(length(e) == size(i, 1) == size(i, 2))
      error("Argument dimensions mismatched")
    end
    return new(e, i)
  end
end
function TransmissionNetworkDistribution(x::Vector{TransmissionNetwork})
  return @inbounds TransmissionNetworkDistribution(mean([y.external for y in x]), mean([y.internal for y in x]))
end

function TransmissionNetworkDistribution(x::Vector{TransmissionNetworkDistribution})
  return @inbounds TransmissionNetworkDistribution(mean([y.external for y in x]), mean([y.internal for y in x]))
end

const TNDistribution = TransmissionNetworkDistribution
const TransmissionNetworkPrior = TransmissionNetworkDistribution
const TNPrior = TransmissionNetworkPrior
const TransmissionNetworkPosterior = TransmissionNetworkDistribution
const TNPosterior = TransmissionNetworkPosterior