struct TransmissionNetworkPosterior
  external::Array{Float64, 1}
  internal::Array{Float64, 2}

  function TransmissionNetworkPosterior(x::Vector{TransmissionNetwork})
    return new(mean([y.external for y in x]), mean([y.internal for y in x]))
  end

  function TransmissionNetworkPosterior(x::Vector{TransmissionNetworkPosterior})
    return new(mean([y.external for y in x]), mean([y.internal for y in x]))
  end
end

const TNPosterior = TransmissionNetworkPosterior

function TransmissionNetworkPosterior(x::MarkovChain)
  return TNPosterior(x.transmission_network)
end

function TransmissionNetworkPosterior(x::MCMC)
  return TNPosterior([TNPosterior(y) for y in x.markov_chains])
end

function Base.sum(x::TNPosterior)
  return sum(x.external) + sum(x.internal)
end

function Base.show(io::IO, x::TNPosterior)
  return print(io, "Transmission network posterior distribution")
end
