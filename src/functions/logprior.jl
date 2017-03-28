"""
logprior(riskpriors::RiskParameterPriors,
         riskparams::RiskParameters)

Calculate the log prior of a set of `RiskParameters`
"""
function logprior(riskpriors::RiskParameterPriors,
                  riskparams::RiskParameters)
  lp = 0.
  for i = 1:length(riskparams.sparks)
    lp += loglikelihood(riskpriors.sparks[i], [riskparams.sparks[i]])
  end
  for i = 1:length(riskparams.susceptibility)
    lp += loglikelihood(riskpriors.susceptibility[i], [riskparams.susceptibility[i]])
  end
  for i = 1:length(riskparams.transmissibility)
    lp += loglikelihood(riskpriors.transmissibility[i], [riskparams.transmissibility[i]])
  end
  for i = 1:length(riskparams.infectivity)
    lp += loglikelihood(riskpriors.infectivity[i], [riskparams.infectivity[i]])
  end
  for i = 1:length(riskparams.latency)
    lp += loglikelihood(riskpriors.latency[i], [riskparams.latency[i]])
  end
  for i = 1:length(riskparams.removal)
    lp += loglikelihood(riskpriors.removal[i], [riskparams.removal[i]])
  end
  return lp
end
