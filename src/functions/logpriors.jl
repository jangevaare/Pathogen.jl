function logpriors(rparams::RiskParameters{T}, rpriors::RiskPriors{T}) where T <: EpidemicModel
  lprior = 0.0
  for i in 1:length(rpriors.sparks)
    lprior += logpdf(rpriors.sparks[i], rparams.sparks[i])
  end
  for i in 1:length(rpriors.susceptibility)
    lprior += logpdf(rpriors.susceptibility[i], rparams.susceptibility[i])
  end
  for i in 1:length(rpriors.transmissibility)
    lprior += logpdf(rpriors.transmissibility[i], rparams.transmissibility[i])
  end
  for i in 1:length(rpriors.infectivity)
    lprior += logpdf(rpriors.infectivity[i], rparams.infectivity[i])
  end
  if T in [SEIR; SEI]
    for i in 1:length(rpriors.latency)
      lprior += logpdf(rpriors.latency[i], rparams.latency[i])
    end
  end
  if T in [SEIR; SIR]
    for i in 1:length(rpriors.removal)
      lprior += logpdf(rpriors.removal[i], rparams.removal[i])
    end
  end
  return lprior
end
