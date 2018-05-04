struct RiskFunctions{T <: EpidemicModel}
  sparks::Function
  susceptibility::Function
  transmissibility::Function
  infectivity::Function
  latency::Function
  removal::Function

  function RiskFunctions{T}(f...)
    return _init_RiskFunctions!(new{T}(), f...)
  end
end

function _init_RiskFunctions!(x::RiskFunctions{SEIR}, f...)
  x.sparks = f[1]
  x.susceptibility = f[2]
  x.transmissibility = f[3]
  x.infectivity = f[4]
  x.latency = f[5]
  x.removal = f[6]
  return x
end

function _init_RiskFunctions!(x::RiskFunctions{SEI}, f...)
  x.sparks = f[1]
  x.susceptibility = f[2]
  x.transmissibility = f[3]
  x.infectivity = f[4]
  x.latency = f[5]
  return x
end

function _init_RiskFunctions!(x::RiskFunctions{SIR}, f...)
  x.sparks = f[1]
  x.susceptibility = f[2]
  x.transmissibility = f[3]
  x.infectivity = f[4]
  x.removal = f[5]
  return x
end

function _init_RiskFunctions!(x::RiskFunctions{SI}, f...)
  x.sparks = f[1]
  x.susceptibility = f[2]
  x.transmissibility = f[3]
  x.infectivity = f[4]
  return x
end
