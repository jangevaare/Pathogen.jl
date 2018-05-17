mutable struct TransmissionRates
  external::Array{Float64, 1}
  internal::Array{Float64, 2}
  # TODO: Add reference back to `EventRates`?

  function TransmissionRates(individuals::Int64)
    return new(fill(0., individuals),
               fill(0., (individuals, individuals)))
  end

  function TransmissionRates(external::Array{Float64, 1},
                             internal::Array{Float64, 2})
    if !(length(external) == size(internal, 1) == size(internal, 2))
      ArgumentError("Argument dimensions mismatched")
    end
    return new(external, internal)
  end
end

function copy(x::TransmissionRates)
  return TransmissionRates(copy(x.external),
                           copy(x.internal))
end
