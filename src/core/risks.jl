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


function length(x::RiskParameters)
  return sum([length(x.sparks);
              length(x.susceptibility);
              length(x.transmissibility);
              length(x.infectivity);
              length(x.latency);
              length(x.removal)])
end


function getindex(x::RiskParameters, i::Real)
  indices = cumsum([length(x.sparks);
                    length(x.susceptibility);
                    length(x.transmissibility);
                    length(x.infectivity);
                    length(x.latency);
                    length(x.removal)])
  riskfunc = findfirst(i <= indices)
  if riskfunc == 0
    throw(BoundsError)
  end
  return x.(fieldnames(x)[riskfunc])[end - (indices[riskfunc] - i)]
end


function getindex(x::Array{RiskParameters, 1}, i::Real, j::Real)
  return x[i][j]
end


function convert(::Type{Array}, x::Array{RiskParameters, 1})
  parameters = length(x[1])
  iterations = length(x)
  return [x[i, j] for i = 1:iterations, j = 1:parameters]
end
