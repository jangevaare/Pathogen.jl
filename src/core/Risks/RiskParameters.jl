struct RiskParameters{S} <: AbstractRisk{S}
  sparks::Union{Nothing, Vector{Float64}}
  susceptibility::Union{Nothing, Vector{Float64}}
  infectivity::Union{Nothing, Vector{Float64}}
  transmissibility::Union{Nothing, Vector{Float64}}
  latency::Union{Nothing, Vector{Float64}}
  removal::Union{Nothing, Vector{Float64}}

  RiskParameters{SEIR}(ϵ, θs, κ, θt, θl, θr) = new{SEIR}(ϵ, θs, κ, θt, θl, θr)
  RiskParameters{SEI}(ϵ, θs, κ, θt, θl) = new{SEI}(ϵ, θs, κ, θt, θl, nothing)
  RiskParameters{SIR}(ϵ, θs, κ, θt, θr) = new{SIR}(ϵ, θs, κ, θt, nothing, θr)
  RiskParameters{SI}(ϵ, θs, κ, θt) = new{SI}(ϵ, θs, κ, θt, nothing, nothing)
end


function Base.show(io::IO, x::RiskParameters{S}) where {S <: DiseaseStateSequence}
  return print(io, "$S model risk function parameters")
end


function Base.convert(::Type{Vector{Float64}},
                      x::RiskParameters)
  return [x[i] for i = 1:length(x)]
end


function Base.convert(::Type{Array{Float64, 2}},
                      x::Vector{RiskParameters{S}}) where {
                      S <: DiseaseStateSequence}
  return [x[i][j] for i = 1:length(x), j = 1:length(x[1])]
end


function Base.similar(x::RiskParameters{S}, 
                      v::Vector{Float64}) where {
                      S <: DiseaseStateSequence}
  indices = _indices(x, zeros=false)
  if indices[end] != length(v)
    throw(ErrorException("Incompatiable parameter vector"))
  end
  if S == SEIR
    return RiskParameters{S}(v[1:(indices[1])],
                             v[(indices[1]+1):(indices[2])],
                             v[(indices[2]+1):(indices[3])],
                             v[(indices[3]+1):(indices[4])],
                             v[(indices[4]+1):(indices[5])],
                             v[(indices[5]+1):(indices[6])])
  elseif S in [SEI; SIR]
    return RiskParameters{S}(v[1:(indices[1])],
                             v[(indices[1]+1):(indices[2])],
                             v[(indices[2]+1):(indices[3])],
                             v[(indices[3]+1):(indices[4])],
                             v[(indices[4]+1):(indices[5])])
  elseif S == SI
    return RiskParameters{S}(v[1:(indices[1])],
                             v[(indices[1]+1):(indices[2])],
                             v[(indices[2]+1):(indices[3])],
                             v[(indices[3]+1):(indices[4])])
  end
end
