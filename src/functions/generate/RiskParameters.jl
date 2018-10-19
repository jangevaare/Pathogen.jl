function generate(::Type{RiskParameters{T}}, rpriors::RiskPriors{T}) where T <: EpidemicModel
  sparks = Float64[rand(x) for x in rpriors.sparks]
  susceptibility = Float64[rand(x) for x in rpriors.susceptibility]
  transmissibility = Float64[rand(x) for x in rpriors.transmissibility]
  infectivity = Float64[rand(x) for x in rpriors.infectivity]
  if T in [SEIR; SEI]
    latency = Float64[rand(x) for x in rpriors.latency]
  end
  if T in [SEIR; SIR]
    removal = Float64[rand(x) for x in rpriors.removal]
  end
  if T == SEIR
    return RiskParameters{T}(sparks,
                             susceptibility,
                             transmissibility,
                             infectivity,
                             latency,
                             removal)
  elseif T == SEI
    return RiskParameters{T}(sparks,
                             susceptibility,
                             transmissibility,
                             infectivity,
                             latency)
  elseif T == SIR
    return RiskParameters{T}(sparks,
                             susceptibility,
                             transmissibility,
                             infectivity,
                             removal)
  elseif T == SI
    return RiskParameters{T}(sparks,
                             susceptibility,
                             transmissibility,
                             infectivity)
  end
end

function generate(::Type{RiskParameters}, mcmc::MCMC{T}) where T <: EpidemicModel
  return generate(RiskParameters{T}, mcmc.risk_priors)
end

function generate(::Type{RiskParameters{T}},
                  last_rparams::RiskParameters{T},
                  Σ::Array{Float64, 2}) where T <: EpidemicModel
  rparams_vector = rand(MvNormal(convert(Vector, last_rparams), Σ))
  return _like(last_rparams, rparams_vector)
end

function generate(::Type{RiskParameters{T}},
                  mc::MarkovChain{T},
                  Σ::Array{Float64, 2}) where T <: EpidemicModel
  return generate(RiskParameters{T}, mc.risk_parameters[end], Σ)
end
