abstract RiskFunctions
abstract RiskParameters


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
Parameter vectors for `SEIR_RiskFunctions`
"""
type SEIR_RiskParameters <: RiskParameters
  sparks::Vector{Float64}
  susceptibility::Vector{Float64}
  transmissibility::Vector{Float64}
  infectivity::Vector{Float64}
  latency::Vector{Float64}
  removal::Vector{Float64}
end


function copy(x::SEIR_RiskParameters)
  return SEIR_RiskParameters(copy(x.sparks),
                             copy(x.susceptibility),
                             copy(x.transmissibility),
                             copy(x.infectivity),
                             copy(x.latency),
                             copy(x.removal))
end


function length(x::SEIR_RiskParameters)
  return sum([length(x.sparks);
              length(x.susceptibility);
              length(x.transmissibility);
              length(x.infectivity);
              length(x.latency);
              length(x.removal)])
end


function getindex(x::SEIR_RiskParameters,
                  i::Int64)
  indices = cumsum([length(x.sparks);
                    length(x.susceptibility);
                    length(x.transmissibility);
                    length(x.infectivity);
                    length(x.latency);
                    length(x.removal)])
  riskfunc = findfirst(i .<= indices)
  return getfield(x, riskfunc)[end - (indices[riskfunc] - i)]
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
Parameter vectors for `SIR_RiskFunctions`
"""
type SIR_RiskParameters <: RiskParameters
  sparks::Vector{Float64}
  susceptibility::Vector{Float64}
  transmissibility::Vector{Float64}
  infectivity::Vector{Float64}
  removal::Vector{Float64}
end


function copy(x::SIR_RiskParameters)
  return SIR_RiskParameters(copy(x.sparks),
                            copy(x.susceptibility),
                            copy(x.transmissibility),
                            copy(x.infectivity),
                            copy(x.removal))
end


function length(x::SIR_RiskParameters)
  return sum([length(x.sparks);
              length(x.susceptibility);
              length(x.transmissibility);
              length(x.infectivity);
              length(x.removal)])
end


function getindex(x::SIR_RiskParameters,
                  i::Int64)
  indices = cumsum([length(x.sparks);
                    length(x.susceptibility);
                    length(x.transmissibility);
                    length(x.infectivity);
                    length(x.removal)])
  riskfunc = findfirst(i .<= indices)
  return getfield(x, riskfunc)[end - (indices[riskfunc] - i)]
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
Parameter vectors for `SEI_RiskFunctions`
"""
type SEI_RiskParameters <: RiskParameters
  sparks::Vector{Float64}
  susceptibility::Vector{Float64}
  transmissibility::Vector{Float64}
  infectivity::Vector{Float64}
  latency::Vector{Float64}
end


function copy(x::SEI_RiskParameters)
  return SEI_RiskParameters(copy(x.sparks),
                            copy(x.susceptibility),
                            copy(x.transmissibility),
                            copy(x.infectivity),
                            copy(x.latency))
end


function length(x::SEI_RiskParameters)
  return sum([length(x.sparks);
              length(x.susceptibility);
              length(x.transmissibility);
              length(x.infectivity);
              length(x.latency)])
end


function getindex(x::SEI_RiskParameters,
                  i::Int64)
  indices = cumsum([length(x.sparks);
                    length(x.susceptibility);
                    length(x.transmissibility);
                    length(x.infectivity);
                    length(x.latency)])
  riskfunc = findfirst(i .<= indices)
  return getfield(x, riskfunc)[end - (indices[riskfunc] - i)]
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


"""
Parameter vectors for `SI_RiskFunctions`
"""
type SI_RiskParameters <: RiskParameters
  sparks::Vector{Float64}
  susceptibility::Vector{Float64}
  transmissibility::Vector{Float64}
  infectivity::Vector{Float64}
end


function copy(x::SI_RiskParameters)
  return SI_RiskParameters(copy(x.sparks),
                           copy(x.susceptibility),
                           copy(x.transmissibility),
                           copy(x.infectivity))
end


function length(x::SI_RiskParameters)
  return sum([length(x.sparks);
              length(x.susceptibility);
              length(x.transmissibility);
              length(x.infectivity)])
end


function getindex(x::SI_RiskParameters,
                  i::Int64)
  indices = cumsum([length(x.sparks);
                    length(x.susceptibility);
                    length(x.transmissibility);
                    length(x.infectivity)])
  riskfunc = findfirst(i .<= indices)
  return getfield(x, riskfunc)[end - (indices[riskfunc] - i)]
end


function convert(::Type{Vector}, x::RiskParameters)
  parameters = length(x)
  return [x[i] for i = 1:parameters]
end


function convert(::Type{Array}, x::Vector{RiskParameters})
  parameters = length(x[1])
  iterations = length(x)
  return [x[i][j] for i = 1:iterations, j = 1:parameters]
end
