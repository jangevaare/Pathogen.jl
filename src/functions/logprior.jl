"""
logprior(riskpriors::SEIR_RiskParameterPriors,
         riskparams::SEIR_RiskParameters)

Calculate the log prior of a set of `SEIR_RiskParameters`
"""
function logprior(riskpriors::SEIR_RiskParameterPriors,
                  riskparams::SEIR_RiskParameters)
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


"""
logprior(riskpriors::SIR_RiskParameterPriors,
         riskparams::SIR_RiskParameters)

Calculate the log prior of a set of `SIR_RiskParameters`
"""
function logprior(riskpriors::SIR_RiskParameterPriors,
                  riskparams::SIR_RiskParameters)
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
  for i = 1:length(riskparams.removal)
    lp += loglikelihood(riskpriors.removal[i], [riskparams.removal[i]])
  end
  return lp
end


"""
logprior(riskpriors::SEI_RiskParameterPriors,
         riskparams::SEI_RiskParameters)

Calculate the log prior of a set of `SEI_RiskParameters`
"""
function logprior(riskpriors::SEI_RiskParameterPriors,
                  riskparams::SEI_RiskParameters)
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
  return lp
end


"""
logprior(riskpriors::SI_RiskParameterPriors,
         riskparams::SI_RiskParameters)

Calculate the log prior of a set of `SI_RiskParameters`
"""
function logprior(riskpriors::SI_RiskParameterPriors,
                  riskparams::SI_RiskParameters)
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
  return lp
end
