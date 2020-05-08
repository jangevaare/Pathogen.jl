function logpriors(rparams::RiskParameters{T}, rpriors::RiskPriors{T}) where T <: DiseaseStateSequence
  lprior = 0.0
  for i in 1:length(rpriors.sparks)
    lprior == -Inf && break
    lprior += logpdf(rpriors.sparks[i], rparams.sparks[i])
  end
  for i in 1:length(rpriors.susceptibility)
    lprior == -Inf && break
    lprior += logpdf(rpriors.susceptibility[i], rparams.susceptibility[i])
  end
  for i in 1:length(rpriors.infectivity)
    lprior == -Inf && break
    lprior += logpdf(rpriors.infectivity[i], rparams.infectivity[i])
  end
  for i in 1:length(rpriors.transmissibility)
    lprior == -Inf && break
    lprior += logpdf(rpriors.transmissibility[i], rparams.transmissibility[i])
  end
  if State_E ∈ T
    for i in 1:length(rpriors.latency)
      lprior == -Inf && break
      lprior += logpdf(rpriors.latency[i], rparams.latency[i])
    end
  end
  if State_R ∈ T
    for i in 1:length(rpriors.removal)
      lprior == -Inf && break
      lprior += logpdf(rpriors.removal[i], rparams.removal[i])
    end
  end
  return lprior
end
