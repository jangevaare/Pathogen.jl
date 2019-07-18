struct RiskFunctions{T<: EpidemicModel}
  sparks::Function
  susceptibility::Function
  infectivity::Function
  transmissibility::Function
  latency::Function
  removal::Function
end

# Placeholder risk function
function Ωx()
  return NaN
end

# Outer constructors
RiskFunctions{SEI}(ϵ, Ωs, Ωi, κ, Ωl) = RiskFunctions{SEI}(ϵ, Ωs, Ωi, κ, Ωl, Ωx)
RiskFunctions{SIR}(ϵ, Ωs, Ωi, κ, Ωr) = RiskFunctions{SIR}(ϵ, Ωs, Ωi, κ, Ωx, Ωr)
RiskFunctions{SI}(ϵ, Ωs, Ωi, κ) = RiskFunctions{SI}(ϵ, Ωs, Ωi, κ, Ωx, Ωx)

function Base.show(io::IO, x::RiskFunctions{T}) where T <: EpidemicModel
  return print(io, "$T model risk functions")
end
