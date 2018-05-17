mutable struct RiskParameters{T <: EpidemicModel}
  sparks::Vector{Float64}
  susceptibility::Vector{Float64}
  transmissibility::Vector{Float64}
  infectivity::Vector{Float64}
  latency::Vector{Float64}
  removal::Vector{Float64}

  function RiskParameters{T}(f...) where T <: EpidemicModel
    return _init_RiskParameters!(new{T}(), f)
  end
end

function _init_RiskParameters!(x::RiskParameters{SEIR}, f)
  if length(f) != 6
    error("Incorrect number of risk parameter vectors for SEIR models")
  end
  x.sparks = f[1]
  x.susceptibility = f[2]
  x.transmissibility = f[3]
  x.infectivity = f[4]
  x.latency = f[5]
  x.removal = f[6]
  return x
end

function _init_RiskParameters!(x::RiskParameters{SEI}, f)
  if length(f) != 5
    error("Incorrect number of risk parameter vectors for SEI models")
  end
  x.sparks = f[1]
  x.susceptibility = f[2]
  x.transmissibility = f[3]
  x.infectivity = f[4]
  x.latency = f[5]
  return x
end

function _init_RiskParameters!(x::RiskParameters{SIR}, f)
  if length(f) != 5
    error("Incorrect number of risk parameter vectors for SIR models")
  end
  x.sparks = f[1]
  x.susceptibility = f[2]
  x.transmissibility = f[3]
  x.infectivity = f[4]
  x.removal = f[5]
  return x
end

function _init_RiskParameters!(x::RiskParameters{SI}, f)
  if length(f) != 4
    error("Incorrect number of risk parameter vectors for SI models")
  end
  x.sparks = f[1]
  x.susceptibility = f[2]
  x.transmissibility = f[3]
  x.infectivity = f[4]
  return x
end

function Base.copy(x::RiskParameters{SEIR})
  return RiskParameters{SEIR}(copy(x.sparks),
                              copy(x.susceptibility),
                              copy(x.transmissibility),
                              copy(x.infectivity),
                              copy(x.latency),
                              copy(x.removal))
end

function Base.copy(x::RiskParameters{SEI})
  return RiskParameters{SEI}(copy(x.sparks),
                             copy(x.susceptibility),
                             copy(x.transmissibility),
                             copy(x.infectivity),
                             copy(x.latency))
end

function Base.copy(x::RiskParameters{SIR})
  return RiskParameters{SIR}(copy(x.sparks),
                             copy(x.susceptibility),
                             copy(x.transmissibility),
                             copy(x.infectivity),
                             copy(x.removal))
end

function Base.copy(x::RiskParameters{SI})
  return RiskParameters{SI}(copy(x.sparks),
                            copy(x.susceptibility),
                            copy(x.transmissibility),
                            copy(x.infectivity))
end

function length(x::RiskParameters{T}) where T <: EpidemicModel
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

function Base.getindex(x::RiskParameters{T},
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

function Base.setindex!(x::RiskParameters{T},
                   z::Float64,
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

function Base.convert(::Type{Vector},
                      x::RiskParameters{T}) where T <: EpidemicModel
  return [x[i] for i = 1:length(x)]
end

function Base.convert(::Type{Array{Float64, 2}},
                      x::Vector{RiskParameters{T}}) where T <: EpidemicModel
  return [x[i][j] for i = 1:length(x), j = 1:length(x[1])]
end
