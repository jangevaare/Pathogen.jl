struct RiskParameters{T <: DiseaseStateSequence, A, B, C, D, E, F} <: AbstractRisk{T}
  sparks::A
  susceptibility::B
  infectivity::C
  transmissibility::D
  latency::E
  removal::F
end

RiskParameters{SEIR}(ϵ::AbstractVector, θs::AbstractVector, κ::AbstractVector, θt::AbstractVector, θl::AbstractVector, θr::AbstractVector) =
  RiskParameters{SEIR, typeof(ϵ), typeof(θs), typeof(κ), typeof(θt), typeof(θl), typeof(θr)}(ϵ, θs, κ, θt, θl, θr)

RiskParameters{SEI}(ϵ::AbstractVector, θs::AbstractVector, κ::AbstractVector, θt::AbstractVector, θl::AbstractVector) =
  RiskParameters{SEI, typeof(ϵ), typeof(θs), typeof(κ), typeof(θt), typeof(θl), Nothing}(ϵ, θs, κ, θt, θl, nothing)

RiskParameters{SIR}(ϵ::AbstractVector, θs::AbstractVector, κ::AbstractVector, θt::AbstractVector, θr::AbstractVector) =
  RiskParameters{SIR, typeof(ϵ), typeof(θs), typeof(κ), typeof(θt), Nothing, typeof(θr)}(ϵ, θs, κ, θt, nothing, θr)

RiskParameters{SI}(ϵ::AbstractVector, θs::AbstractVector, κ::AbstractVector, θt::AbstractVector) =
  RiskParameters{SI, typeof(ϵ), typeof(θs), typeof(κ), typeof(θt), Nothing, Nothing}(ϵ, θs, κ, θt, nothing, nothing)

function Base.show(io::IO, x::RiskParameters{T}) where T <: DiseaseStateSequence
  return print(io, "$T model risk function parameters")
end

function Base.convert(::Type{Vector{Float64}},
                      x::RiskParameters{T}) where T <: DiseaseStateSequence
  return [x[i] for i = 1:length(x)]
end

function Base.convert(::Type{Array{Float64, 2}},
                      x::Vector{RiskParameters{T}}) where T <: DiseaseStateSequence
  return [x[i][j] for i = 1:length(x), j = 1:length(x[1])]
end

function _like(x::RiskParameters{T}, v::Vector{Float64}) where T <: DiseaseStateSequence
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
  elseif T in (SEI, SIR)
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
