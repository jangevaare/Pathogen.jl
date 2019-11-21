struct RiskParameters{T} <: AbstractRisk{T}
  sparks::Union{Nothing, AbstractVector}
  susceptibility::Union{Nothing, AbstractVector}
  infectivity::Union{Nothing, AbstractVector}
  transmissibility::Union{Nothing, AbstractVector}
  latency::Union{Nothing, AbstractVector}
  removal::Union{Nothing, AbstractVector}

  RiskParameters{SEIR}(ϵ, θs, θi, κ, θl, θr) = new{SEIR}(ϵ, θs, θi, κ, θl, θr)
  RiskParameters{SEI}(ϵ, θs, θi, κ, θl) = new{SEI}(ϵ, θs, θi, κ, θl, nothing)
  RiskParameters{SIR}(ϵ, θs, θi, κ, θr) = new{SIR}(ϵ, θs, θi, κ, nothing, θr)
  RiskParameters{SI}(ϵ, θs, θi, κ) = new{SI}(ϵ, θs, θi, κ, nothing, nothing)
end


function Base.show(io::IO, x::RiskParameters{T}) where T <: EpidemicModel
  return print(io, "$T model risk function parameters")
end


function Base.convert(::Type{Vector{Float64}},
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
    @error "Incompatiable parameter vector"
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
