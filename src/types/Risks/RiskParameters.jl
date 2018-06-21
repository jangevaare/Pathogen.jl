mutable struct RiskParameters{T <: EpidemicModel}
  sparks::Vector{Float64}
  susceptibility::Vector{Float64}
  transmissibility::Vector{Float64}
  infectivity::Vector{Float64}
  latency::Vector{Float64}
  removal::Vector{Float64}

  function RiskParameters{T}(sp::V, su::V, tr::V, in::V, la::V, re::V) where
    {T <: SEIR, V <: Vector{Float64}}
    return new{T}(sp, su, tr, in, la, re)
  end

  function RiskParameters{T}(sp::V, su::V, tr::V, in::V, la::V) where
    {T <: SEI, V <: Vector{Float64}}
    x = new{T}()
    x.sparks = sp
    x.susceptibility = su
    x.transmissibility = tr
    x.infectivity = in
    x.latency = la
    return x
  end

  function RiskParameters{T}(sp::V, su::V, tr::V, in::V, re::V) where
    {T <: SIR, V <: Vector{Float64}}
    x = new{T}()
    x.sparks = sp
    x.susceptibility = su
    x.transmissibility = tr
    x.infectivity = in
    x.removal = re
    return x
  end

  function RiskParameters{T}(sp::V, su::V, tr::V, in::V) where
    {T <: SI, V <: Vector{Float64}}
    x = new{T}()
    x.sparks = sp
    x.susceptibility = su
    x.transmissibility = tr
    x.infectivity = in
    return x
  end
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

function _indices(x::RiskParameters{T}; zeros::Bool=true) where T <: EpidemicModel
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

function length(x::RiskParameters{T}) where T <: EpidemicModel
  return _indices(x)[end]
end

function Base.getindex(x::RiskParameters{T},
                       i::Int64) where T <: EpidemicModel
  indices = _indices(x, zeros = true)
  riskfunc = findfirst(i .<= indices)
  return getfield(x, riskfunc)[end - (indices[riskfunc] - i)]
end

function Base.setindex!(x::RiskParameters{T},
                        z::Float64,
                        i::Int64) where T <: EpidemicModel
  indices = _indices(x, zeros = true)
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

function _like(x::RiskParameters{T}, v::Vector{Float64}) where T <: EpidemicModel
  indices = _indices(x, zeros=false)
  if indices[end] != length(v)
    error("Incompatiable parameter vector")
  end
  if T == SEIR
    return RiskParameters{T}(v[1:(indices[1])],
                             v[(indices[1]+1):(indices[2])],
                             v[(indices[2]+1):(indices[3])],
                             v[(indices[3]+1):(indices[4])],
                             v[(indices[4]+1):(indices[5])],
                             v[(indices[5]+1):(indices[6])])
  elseif T in [SEI; SIR]
    return RiskParameters{T}(v[1:(indices[1])],
                             v[(indices[1]+1):(indices[2])],
                             v[(indices[2]+1):(indices[3])],
                             v[(indices[3]+1):(indices[4])],
                             v[(indices[4]+1):(indices[5])])
  elseif T == SI
    return RiskParameters{T}(v[1:(indices[1])],
                             v[(indices[1]+1):(indices[2])],
                             v[(indices[2]+1):(indices[3])],
                             v[(indices[3]+1):(indices[4])])
  end
end
