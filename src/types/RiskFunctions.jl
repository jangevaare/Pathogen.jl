abstract RiskFunctions


"""
A collection of the risk functions required for SEIR simulation and inference
"""
type SEIR_RiskFunctions <: RiskFunctions
  sparks::Function
  susceptibility::Function
  transmissibility::Function
  infectivity::Function
  latency::Function
  removal::Function
end


"""
A collection of the risk functions required for SIR simulation and inference
"""
type SIR_RiskFunctions <: RiskFunctions
  sparks::Function
  susceptibility::Function
  transmissibility::Function
  infectivity::Function
  removal::Function
end


"""
A collection of the risk functions required for SEI simulation and inference
"""
type SEI_RiskFunctions <: RiskFunctions
  sparks::Function
  susceptibility::Function
  transmissibility::Function
  infectivity::Function
  latency::Function
end


"""
A collection of the risk functions required for SI simulation and inference
"""
type SI_RiskFunctions <: RiskFunctions
  sparks::Function
  susceptibility::Function
  transmissibility::Function
  infectivity::Function
end
