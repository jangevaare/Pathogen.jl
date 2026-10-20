struct RiskFunctions{T <: DiseaseStateSequence, A, B, C, D, E, F}
  sparks::A
  susceptibility::B
  infectivity::C
  transmissibility::D
  latency::E
  removal::F
end

RiskFunctions{SEIR}(ϵ, Ωs, κ, Ωt, Ωl, Ωr) =
  RiskFunctions{SEIR, typeof(ϵ), typeof(Ωs), typeof(κ), typeof(Ωt), typeof(Ωl), typeof(Ωr)}(ϵ, Ωs, κ, Ωt, Ωl, Ωr)

RiskFunctions{SEI}(ϵ, Ωs, κ, Ωt, Ωl) =
  RiskFunctions{SEI, typeof(ϵ), typeof(Ωs), typeof(κ), typeof(Ωt), typeof(Ωl), Nothing}(ϵ, Ωs, κ, Ωt, Ωl, nothing)

RiskFunctions{SIR}(ϵ, Ωs, κ, Ωt, Ωr) =
  RiskFunctions{SIR, typeof(ϵ), typeof(Ωs), typeof(κ), typeof(Ωt), Nothing, typeof(Ωr)}(ϵ, Ωs, κ, Ωt, nothing, Ωr)

RiskFunctions{SI}(ϵ, Ωs, κ, Ωt) =
  RiskFunctions{SI, typeof(ϵ), typeof(Ωs), typeof(κ), typeof(Ωt), Nothing, Nothing}(ϵ, Ωs, κ, Ωt, nothing, nothing)

function Base.show(io::IO, x::RiskFunctions{T}) where T <: DiseaseStateSequence
  return print(io, "$T model risk functions")
end
