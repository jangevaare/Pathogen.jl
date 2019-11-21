struct RiskFunctions{T} <: AbstractRisk{T}
  sparks::Union{Nothing, Function}
  susceptibility::Union{Nothing, Function}
  infectivity::Union{Nothing, Function}
  transmissibility::Union{Nothing, Function}
  latency::Union{Nothing, Function}
  removal::Union{Nothing, Function}

  RiskFunctions{SEIR}(ϵ, Ωs, Ωi, κ, Ωl, Ωr) = new{SEIR}(ϵ, Ωs, Ωi, κ, Ωl, Ωr)
  RiskFunctions{SEI}(ϵ, Ωs, Ωi, κ, Ωl) = new{SEI}(ϵ, Ωs, Ωi, κ, Ωl, nothing)
  RiskFunctions{SIR}(ϵ, Ωs, Ωi, κ, Ωr) = new{SIR}(ϵ, Ωs, Ωi, κ, nothing, Ωr)
  RiskFunctions{SI}(ϵ, Ωs, Ωi, κ) = new{SI}(ϵ, Ωs, Ωi, κ, nothing, nothing)
end


function Base.show(io::IO, x::RiskFunctions{T}) where T <: EpidemicModel
  return print(io, "$T model risk functions")
end
