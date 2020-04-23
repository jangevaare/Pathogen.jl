struct TransmissionRates <: AbstractTN
  external::Vector{Float64}
  internal::Array{Float64, 2}

  function TransmissionRates(individuals::Integer)
    return new(fill(0., individuals),
               fill(0., (individuals, individuals)))
  end

  function TransmissionRates(e::Vector{Float64},
                             i::Array{Float64, 2})
    @boundscheck if !(length(e) == size(i, 1) == size(i, 2))
      error("Argument dimensions mismatched")
    end
    return new(e, i)
  end
end