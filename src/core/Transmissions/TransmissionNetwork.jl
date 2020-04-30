struct TransmissionNetwork <: AbstractTN
  external::BitArray{1}
  internal::BitArray{2}

  function TransmissionNetwork(individuals::Int64)
    return new(fill(0, individuals),
               fill(0, (individuals, individuals)))
  end

  function TransmissionNetwork(external::BitArray{1},
                               internal::BitArray{2})
    if !(length(external) == size(internal, 1) == size(internal, 2))
      throw(ErrorException("Mismatched BitArray sizes"))
    end
    multiple_exposures = findall((sum(internal, dims=1)[:] .+ external) .> 1)
    if length(multiple_exposures) > 0
      throw(ErrorException("Multiple exposures detected for individual(s): $multiple_exposures"))
    end
    return new(external, internal)
  end

  function TransmissionNetwork(starting_states::DiseaseStates)
    return new(starting_states .!= Ref(State_S), fill(0, (length(starting_states), length(starting_states))))
  end
end