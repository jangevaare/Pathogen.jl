"""
A collection of the risk functions required for simulation and inference. These
risk functions are parameterized by a `RiskParameters` object
"""
type RiskFunctions
  susceptibility::Function
  transmissibility::Function
  infectivity::Function
  sparks::Function
  latency::Function
  removal::Function
  detection::Function
end


"""
Parameter vectors for the `RiskFunctions`
"""
type RiskParameters
  susceptibility::Vector{Float64}
  transmissibility::Vector{Float64}
  infectivity::Vector{Float64}
  sparks::Vector{Float64}
  latency::Vector{Float64}
  removal::Vector{Float64}
  detection::Vector{Float64}
end
