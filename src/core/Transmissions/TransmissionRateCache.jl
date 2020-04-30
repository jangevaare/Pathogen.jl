struct TransmissionRateCache <: AbstractTN
  external::Vector{Union{Nothing, Float64}}
  internal::Array{Union{Nothing, Float64}, 2}

  function TransmissionRateCache(x::Integer)
    return new(Union{Nothing, Float64}[nothing for i = 1:x],
               Union{Nothing, Float64}[nothing for i = 1:x, j = 1:x])
  end
end