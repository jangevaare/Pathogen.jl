abstract type RiskParameters end


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


function length(x::RiskParameters)
  params = sum([length(x.sparks);
                length(x.susceptibility);
                length(x.transmissibility);
                length(x.infectivity)])
  if typeof(x) in [SEIR_RiskParameters;
                   SEI_RiskParameters]
    params += length(x.latency)
  end
  if typeof(x) in [SEIR_RiskParameters;
                   SIR_RiskParameters]
    params += length(x.removal)
  end
  return params
end


function getindex(x::RiskParameters,
                  i::Int64)
  indices = Int64[]
  if typeof(x) == SEIR_RiskParameters
    indices = cumsum([length(x.sparks);
                      length(x.susceptibility);
                      length(x.transmissibility);
                      length(x.infectivity);
                      length(x.latency);
                      length(x.removal)])
  elseif typeof(x) == SIR_RiskParameters
    indices = cumsum([length(x.sparks);
                      length(x.susceptibility);
                      length(x.transmissibility);
                      length(x.infectivity);
                      length(x.removal)])
  elseif typeof(x) == SEI_RiskParameters
    indices = cumsum([length(x.sparks);
                      length(x.susceptibility);
                      length(x.transmissibility);
                      length(x.infectivity);
                      length(x.latency)])
  elseif typeof(x) == SI_RiskParameters
    indices = cumsum([length(x.sparks);
                      length(x.susceptibility);
                      length(x.transmissibility);
                      length(x.infectivity)])
  else
    error("Unrecognized Risk Parameter type")
  end
  riskfunc = findfirst(i .<= indices)
  return getfield(x, riskfunc)[end - (indices[riskfunc] - i)]
end


function setindex!(x::RiskParameters,
                   z::Float64,
                   i::Int64)
  indices = Int64[]
  if typeof(x) == SEIR_RiskParameters
    indices = cumsum([length(x.sparks);
                      length(x.susceptibility);
                      length(x.transmissibility);
                      length(x.infectivity);
                      length(x.latency);
                      length(x.removal)])
  elseif typeof(x) == SIR_RiskParameters
    indices = cumsum([length(x.sparks);
                      length(x.susceptibility);
                      length(x.transmissibility);
                      length(x.infectivity);
                      length(x.removal)])
  elseif typeof(x) == SEI_RiskParameters
    indices = cumsum([length(x.sparks);
                      length(x.susceptibility);
                      length(x.transmissibility);
                      length(x.infectivity);
                      length(x.latency)])
  elseif typeof(x) == SI_RiskParameters
    indices = cumsum([length(x.sparks);
                      length(x.susceptibility);
                      length(x.transmissibility);
                      length(x.infectivity)])
  else
    error("Unrecognized Risk Parameter type")
  end
  riskfunc = findfirst(i .<= indices)
  getfield(x, riskfunc)[end - (indices[riskfunc] - i)] = z
  return x
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
