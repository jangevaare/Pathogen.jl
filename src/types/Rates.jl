"""
Network disease transmission rates
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


function copy(x::NetworkRates)
  return NetworkRates(copy(x.external),
                      copy(x.internal))
end


abstract Rates


"""
SEIR event rates
"""
type SEIR_Rates <: Rates
  exposure::NetworkRates
  infection::Vector{Float64}
  removal::Vector{Float64}
  individuals::Int64

  function SEIR_Rates(individuals::Int64)
    return new(NetworkRates(individuals),
               fill(0., individuals),
               fill(0., individuals),
               individuals)
  end
end


function getindex(x::SEIR_Rates, i::Int64)
  if i == 1
    return x.exposure
  elseif i == 2
    return x.infection
  elseif i == 3
    return x.removal
  else
    throw(BoundsError)
  end
end


"""
SIR event rates
"""
type SIR_Rates <: Rates
  infection::NetworkRates
  removal::Vector{Float64}
  individuals::Int64

  function SIR_Rates(individuals::Int64)
    return new(NetworkRates(individuals),
               fill(0., individuals),
               individuals)
  end
end


function getindex(x::SIR_Rates, i::Int64)
  if i == 1
    return x.infection
  elseif i == 2
    return x.removal
  else
    throw(BoundsError)
  end
end


"""
SEI event rates
"""
type SEI_Rates <: Rates
  exposure::NetworkRates
  infection::Vector{Float64}
  individuals::Int64

  function SEI_Rates(individuals::Int64)
    return new(NetworkRates(individuals),
               fill(0., individuals),
               individuals)
  end
end


function getindex(x::SEI_Rates, i::Int64)
  if i == 1
    return x.exposure
  elseif i == 2
    return x.infection
  else
    throw(BoundsError)
  end
end


"""
SI event rates
"""
type SI_Rates <: Rates
  infection::NetworkRates
  individuals::Int64

  function SI_Rates(individuals::Int64)
    return new(NetworkRates(individuals),
               individuals)
  end
end


function getindex(x::SI_Rates, i::Int64)
  if i == 1
    return x.infection
  else
    throw(BoundsError)
  end
end
