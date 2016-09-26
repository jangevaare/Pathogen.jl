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


function size(x::Vector{RiskParameters})
  return (length(x), length(x[1]))
end


function getindex(x::RiskParameters, i)
  inds = cumsum([length(x.sparks);
                 length(x.susceptibility);
                 length(x.transmissibility);
                 length(x.infectivity);
                 length(x.latency);
                 length(x.removal)])
  riskfunc = findfirst(i <= inds)
  if riskfunc == 0
    @error("BoundsError")
  end
  return x.(fieldnames(x)[riskfunc])[end - (inds[riskfunc] - i)]
end


function Vector(x::RiskParameters)
  return [x[i] for i = 1:length(x)]
end


function getindex(x::Vector{RiskParameters}, i, j)
  inds = cumsum([length(x[i].sparks);
                 length(x[i].susceptibility);
                 length(x[i].transmissibility);
                 length(x[i].infectivity);
                 length(x[i].latency);
                 length(x[i].removal)])
  riskfunc = findfirst(j <= inds)
  if riskfunc == 0
    @error("BoundsError")
  end
  return x[i].(fieldnames(x[i])[riskfunc])[end - (inds[riskfunc] - j)]
end


function Array(x::Vector{RiskParameters})
  dims = size(x)
  return [x[i, j] for i = 1:dim[1], j = 1:dim[2]]
end
