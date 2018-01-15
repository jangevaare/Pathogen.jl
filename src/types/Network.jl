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
    multiple_exposures = find((sum(internal, 1)[:] .+ external) .> 1)
    if length(multiple_exposures) > 0
      error("Multiple exposures detected for individual(s): $multiple_exposures")
    end
    return new(external, internal)
  end
end


function copy(x::Network)
  return Network(copy(x.external),
                 copy(x.internal))
end


function show(io::IO, object::Network)
  print(io, "Transmission network with $(sum(object.external)) external
  exposure(s), and $(sum(object.internal)) internal exposure(s)")
end
