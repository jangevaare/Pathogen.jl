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
    elseif any([any(internal[:,i]) for i = 1:length(external)] .& external)
      error("Multiple exposures per individual detected")
    end
    return new(external, internal)
  end
end


function copy(x::Network)
  return Network(copy(x.external),
                 copy(x.internal))
end
