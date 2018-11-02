mutable struct TransmissionNetwork
  external::Vector{Bool}
  internal::Array{Bool, 2}

  function TransmissionNetwork(individuals::Int64)
    return new(fill(false, individuals),
               fill(false, (individuals, individuals)))
  end

  function TransmissionNetwork(external::Vector{Bool},
                               internal::Array{Bool, 2})
    if !(length(external) == size(internal, 1) == size(internal, 2))
      throw(BoundsError)
    end
    multiple_exposures = findall((sum(internal, 1)[:] .+ external) .> 1)
    if length(multiple_exposures) > 0
      @error "Multiple exposures detected for individual(s): $multiple_exposures"
    end
    return new(external, internal)
  end
end

function Base.copy(x::TransmissionNetwork)
  return TransmissionNetwork(copy(x.external),
                             copy(x.internal))
end

function Base.show(io::IO, object::TransmissionNetwork)
  return print(io, "Transmission network with $(sum(object.external)) external, and $(sum(object.internal)) internal transmission(s)")
end
