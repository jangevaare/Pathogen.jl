mutable struct RiskPriors{T <: EpidemicModel}
  sparks::Vector{UnivariateDistribution}
  susceptibility::Vector{UnivariateDistribution}
  transmissibility::Vector{UnivariateDistribution}
  infectivity::Vector{UnivariateDistribution}
  latency::Vector{UnivariateDistribution}
  removal::Vector{UnivariateDistribution}

  function RiskPriors{T}(sp, su, tr, in, la, re) where {T <: SEIR}
    return new{T}(sp, su, tr, in, la, re)
  end

  function RiskPriors{T}(sp, su, tr, in, la) where {T <: SEI}
    x = new{T}()
    x.sparks = sp
    x.susceptibility = su
    x.transmissibility = tr
    x.infectivity = in
    x.latency = la
    return x
  end

  function RiskPriors{T}(sp, su, tr, in, re) where {T <: SIR}
    x = new{T}()
    x.sparks = sp
    x.susceptibility = su
    x.transmissibility = tr
    x.infectivity = in
    x.removal = re
    return x
  end

  function RiskPriors{T}(sp, su, tr, in) where {T <: SI}
    x = new{T}()
    x.sparks = sp
    x.susceptibility = su
    x.transmissibility = tr
    x.infectivity = in
    return x
  end
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

function _indices(x::RiskPriors{T}; zeros::Bool=true) where T <: EpidemicModel
  indices = [length(x.sparks)
             length(x.susceptibility)
             length(x.transmissibility)
             length(x.infectivity)]
  if T in [SEIR; SEI]
    push!(indices, length(x.latency))
  elseif zeros
    push!(indices, 0)
  end
  if T in [SEIR; SIR]
    push!(indices, length(x.removal))
  elseif zeros
    push!(indices, 0)
  end
  return cumsum(indices)
end

function length(x::RiskPriors{T}) where T <: EpidemicModel
  return _indices(x)[end]
end

function Base.getindex(x::RiskPriors{T},
                       i::Int64) where T <: EpidemicModel
  indices = _indices(x, zeros = true)
  riskfunc = findfirst(i .<= indices)
  return getfield(x, riskfunc)[end - (indices[riskfunc] - i)]
end

function Base.setindex!(x::RiskPriors{T},
                        z::UnivariateDistribution,
                        i::Int64) where T <: EpidemicModel
  indices = _indices(x, zeros = true)
  riskfunc = findfirst(i .<= indices)
  getfield(x, riskfunc)[end - (indices[riskfunc] - i)] = z
  return x
end
