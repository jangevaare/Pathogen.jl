struct RiskFunctions{T <: EpidemicModel}
  sparks::Function
  susceptibility::Function
  transmissibility::Function
  infectivity::Function
  latency::Function
  removal::Function

  function RiskFunctions{T}(f...) where T <: EpidemicModel
    return _init_RiskFunctions!(new{T}(), f)
  end
end

function _init_RiskFunctions!(x::RiskFunctions{SEIR}, f)
  if length(f) != 6
    error("Incorrect number of risk functions for SEIR models")
  end
  x.sparks = f[1]
  x.susceptibility = f[2]
  x.transmissibility = f[3]
  x.infectivity = f[4]
  x.latency = f[5]
  x.removal = f[6]
  return x
end

function _init_RiskFunctions!(x::RiskFunctions{SEI}, f)
  if length(f) != 5
    error("Incorrect number of risk functions for SEI models")
  end
  x.sparks = f[1]
  x.susceptibility = f[2]
  x.transmissibility = f[3]
  x.infectivity = f[4]
  x.latency = f[5]
  return x
end

function _init_RiskFunctions!(x::RiskFunctions{SIR}, f)
  if length(f) != 5
    error("Incorrect number of risk functions for SIR models")
  end
  x.sparks = f[1]
  x.susceptibility = f[2]
  x.transmissibility = f[3]
  x.infectivity = f[4]
  x.removal = f[5]
  return x
end

function _init_RiskFunctions!(x::RiskFunctions{SI}, f)
  if length(f) != 6
    error("Incorrect number of risk functions for SI models")
  end
  x.sparks = f[1]
  x.susceptibility = f[2]
  x.transmissibility = f[3]
  x.infectivity = f[4]
  return x
end
