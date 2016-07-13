"""
A collection of the risk functions required for simulation and inference. These
risk functions are parameterized by a `RiskParameters` object
"""
type RiskFunctions
  sparks::Function
  susceptibility::Function
  transmissibility::Function
  infectivity::Function
  latency::Function
  removal::Function
end


"""
Parameter vectors for the `RiskFunctions`
"""
type RiskParameters
  sparks::Vector{Float64}
  susceptibility::Vector{Float64}
  transmissibility::Vector{Float64}
  infectivity::Vector{Float64}
  latency::Vector{Float64}
  removal::Vector{Float64}
end
