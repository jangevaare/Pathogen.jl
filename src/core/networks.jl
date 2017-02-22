"""
An infectious disease transmission network
"""
type Network
  external::Vector{Bool}
  internal::Array{Bool, 2}

  function Network(population::DataFrame)
    individuals = size(population, 1)
    return new(fill(false, individuals),
               fill(false, (individuals, individuals)))
  end

  function Network(individuals::Int64)
    return new(fill(false, individuals),
               fill(false, (individuals, individuals)))
  end

  function Network(external::Vector{Bool},
                   internal::Array{Bool, 2})
    if !(length(external) == size(internal, 1) == size(internal, 2))
      throw(BoundsError)
    end
    return new(external, internal)
  end
end


function copy(network::Network)
  return Network(copy(network.external),
                 copy(network.internal))
end


"""
Exposure rates
"""
type NetworkRates
  external::Vector{Float64}
  internal::Array{Float64, 2}

  function NetworkRates(individuals::Int64)
    return new(fill(0., individuals),
               fill(0., (individuals, individuals)))
  end

  function NetworkRates(external::Vector{Float64},
                        internal::Array{Float64, 2})
    if !(length(external) == size(internal, 1) == size(internal, 2))
      throw(BoundsError)
    end
    return new(external, internal)
  end
end
