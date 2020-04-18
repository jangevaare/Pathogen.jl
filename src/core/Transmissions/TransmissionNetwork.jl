struct TransmissionNetwork
  external::BitArray{1}
  internal::BitArray{2}

  function TransmissionNetwork(individuals::Integer)
    return new(fill(0, individuals),
               fill(0, (individuals, individuals)))
  end

  function TransmissionNetwork(external::BitArray{1},
                               internal::BitArray{2})
    if !(length(external) == size(internal, 1) == size(internal, 2))
      @error "Mismatched BitArray sizes"
    end
    multiple_exposures = findall((sum(internal, dims=1)[:] .+ external) .> 1)
    if length(multiple_exposures) > 0
      @error "Multiple exposures detected for individual(s): $multiple_exposures"
    end
    return new(external, internal)
  end

  function TransmissionNetwork(starting_states::Vector{DiseaseState})
    return new(starting_states .!= Ref(State_S), fill(0, (length(starting_states), length(starting_states))))
  end
end

function Base.copy(x::TransmissionNetwork)
  return TransmissionNetwork(copy(x.external),
                             copy(x.internal))
end

function Base.show(io::IO, object::TransmissionNetwork)
  return print(io, "Transmission network with $(sum(object.external)) external, and $(sum(object.internal)) internal transmission(s)")
end
