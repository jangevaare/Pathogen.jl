"""
Event rates
"""
type Rates
  external::Vector{Float64}
  internal::Array{Float64, 2}
  infection::Vector{Float64}
  removal::Vector{Float64}
  individuals::Int64

  function Rates(individuals::Int64)
    return new(fill(0., individuals),
               fill(0., (individuals, individuals)),
               fill(0., individuals),
               fill(0., individuals),
               individuals)
  end
end


function getindex(x::Rates, i::Int64)
  if i == 1
    return x.external
  elseif i == 2
    return x.internal
  elseif i == 3
    return x.infection
  elseif i == 4
    return x.removal
  else
    throw(BoundsError)
  end
end
