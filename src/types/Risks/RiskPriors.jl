struct RiskPriors{T} where T <: EpidemicModel
  sparks::Vector{UnivariateDistribution}
  susceptibility::Vector{UnivariateDistribution}
  transmissibility::Vector{UnivariateDistribution}
  infectivity::Vector{UnivariateDistribution}
  latency::Vector{UnivariateDistribution}
  removal::Vector{UnivariateDistribution}

  function RiskPriors{T}(f...)
    return _init_RiskPriors!(new{T}(), f...)
  end
end

function _init_RiskPriors!(x::RiskPriors{SEIR}, f...)
  x.sparks = f[1]
  x.susceptibility = f[2]
  x.transmissibility = f[3]
  x.infectivity = f[4]
  x.latency = f[5]
  x.removal = f[6]
  return x
end

function _init_RiskPriors!(x::RiskPriors{SEI}, f...)
  x.sparks = f[1]
  x.susceptibility = f[2]
  x.transmissibility = f[3]
  x.infectivity = f[4]
  x.latency = f[5]
  return x
end

function _init_RiskPriors!(x::RiskPriors{SIR}, f...)
  x.sparks = f[1]
  x.susceptibility = f[2]
  x.transmissibility = f[3]
  x.infectivity = f[4]
  x.removal = f[5]
  return x
end

function _init_RiskPriors!(x::RiskPriors{SI}, f...)
  x.sparks = f[1]
  x.susceptibility = f[2]
  x.transmissibility = f[3]
  x.infectivity = f[4]
  return x
end

function Base.copy(x::RiskPriors{SEIR})
  return RiskPriors{SEIR}(copy(x.sparks),
                              copy(x.susceptibility),
                              copy(x.transmissibility),
                              copy(x.infectivity),
                              copy(x.latency),
                              copy(x.removal))
end

function Base.copy(x::RiskPriors{SEI})
  return RiskPriors{SEI}(copy(x.sparks),
                             copy(x.susceptibility),
                             copy(x.transmissibility),
                             copy(x.infectivity),
                             copy(x.latency))
end

function Base.copy(x::RiskPriors{SIR})
  return RiskPriors{SIR}(copy(x.sparks),
                             copy(x.susceptibility),
                             copy(x.transmissibility),
                             copy(x.infectivity),
                             copy(x.removal))
end

function Base.copy(x::RiskPriors{SI})
  return RiskPriors{SI}(copy(x.sparks),
                            copy(x.susceptibility),
                            copy(x.transmissibility),
                            copy(x.infectivity))
end

function length(x::RiskPriors{T}) where T <: EpidemicModel
  params = sum([length(x.sparks)
                length(x.susceptibility)
                length(x.transmissibility)
                length(x.infectivity)])
  if T in [SEIR; SEI]
    params += length(x.latency)
  end
  if T in [SEIR; SIR]
    params += length(x.removal)
  end
  return params
end

function Base.getindex(x::RiskPriors{T},
                  i::Int64) where T <: EpidemicModel
  indices = [length(x.sparks)
             length(x.susceptibility)
             length(x.transmissibility)
             length(x.infectivity)]
  if T in [SEIR; SEI]
    push!(indices, length(x.latency))
  end
  if T in [SEIR; SIR]
    push!(indices, length(x.removal))
  end
  indices = cumsum(indices)
  riskfunc = findfirst(i .<= indices)
  return getfield(x, riskfunc)[end - (indices[riskfunc] - i)]
end

function Base.setindex!(x::RiskPriors{T},
                   z::UnivariateDistribution,
                   i::Int64) where T <: EpidemicModel
  indices = [length(x.sparks)
             length(x.susceptibility)
             length(x.transmissibility)
             length(x.infectivity)]
  if T in [SEIR; SEI]
    push!(indices, length(x.latency))
  end
  if T in [SEIR; SIR]
    push!(indices, length(x.removal))
  end
  indices = cumsum(indices)
  riskfunc = findfirst(i .<= indices)
  getfield(x, riskfunc)[end - (indices[riskfunc] - i)] = z
  return x
end
