
function generate(::Type{RiskParameters{S}}, 
                  rpriors::RiskPriors{S}) where {
                  S <: DiseaseStateSequence}
  sparks = Float64[rand(x) for x in rpriors.sparks]
  susceptibility = Float64[rand(x) for x in rpriors.susceptibility]
  infectivity = Float64[rand(x) for x in rpriors.infectivity]
  transmissibility = Float64[rand(x) for x in rpriors.transmissibility]
  if S in [SEIR; SEI]
    latency = Float64[rand(x) for x in rpriors.latency]
  end
  if S in [SEIR; SIR]
    removal = Float64[rand(x) for x in rpriors.removal]
  end
  if S == SEIR
    return RiskParameters{S}(sparks,
                             susceptibility,
                             infectivity,
                             transmissibility,
                             latency,
                             removal)
  elseif S == SEI
    return RiskParameters{S}(sparks,
                             susceptibility,
                             infectivity,
                             transmissibility,
                             latency)
  elseif S == SIR
    return RiskParameters{S}(sparks,
                             susceptibility,
                             infectivity,
                             transmissibility,
                             removal)
  elseif S == SI
    return RiskParameters{S}(sparks,
                             susceptibility,
                             infectivity,
                             transmissibility)
  end
end

function generate(::Type{RiskParameters}, mcmc::MCMC{S, M}) where {
                  S <: DiseaseStateSequence, 
                  M <: ILM}
  return generate(RiskParameters{S}, mcmc.risk_priors)
end

function generate(::Type{RiskParameters{S}},
                  last_rparams::RiskParameters{S},
                  Σ::Array{Float64, 2}) where {
                  S <: DiseaseStateSequence}
  rparams_vector = rand(MvNormal(convert(Vector{Float64}, last_rparams), Σ))
  return similar(last_rparams, rparams_vector)
end

function generate(::Type{RiskParameters{S}},
                  mc::MarkovChain{S, M},
                  Σ::Array{Float64, 2}) where {
                  S <: DiseaseStateSequence, 
                  M <: ILM}
  return generate(RiskParameters{S}, mc.risk_parameters[end], Σ)
end