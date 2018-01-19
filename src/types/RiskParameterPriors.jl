abstract type RiskParameterPriors end


"""
Prior distributions vectors for the `SEIR_RiskParameters`
"""
type SEIR_RiskParameterPriors <: RiskParameterPriors
  sparks::Vector{UnivariateDistribution}
  susceptibility::Vector{UnivariateDistribution}
  transmissibility::Vector{UnivariateDistribution}
  infectivity::Vector{UnivariateDistribution}
  latency::Vector{UnivariateDistribution}
  removal::Vector{UnivariateDistribution}
end


"""
Prior distributions vectors for the `SIR_RiskParameters`
"""
type SIR_RiskParameterPriors <: RiskParameterPriors
  sparks::Vector{UnivariateDistribution}
  susceptibility::Vector{UnivariateDistribution}
  transmissibility::Vector{UnivariateDistribution}
  infectivity::Vector{UnivariateDistribution}
  removal::Vector{UnivariateDistribution}
end


"""
Prior distributions vectors for the `SEI_RiskParameters`
"""
type SEI_RiskParameterPriors <: RiskParameterPriors
  sparks::Vector{UnivariateDistribution}
  susceptibility::Vector{UnivariateDistribution}
  transmissibility::Vector{UnivariateDistribution}
  infectivity::Vector{UnivariateDistribution}
  latency::Vector{UnivariateDistribution}
end


"""
Prior distributions vectors for the `SI_RiskParameters`
"""
type SI_RiskParameterPriors <: RiskParameterPriors
  sparks::Vector{UnivariateDistribution}
  susceptibility::Vector{UnivariateDistribution}
  transmissibility::Vector{UnivariateDistribution}
  infectivity::Vector{UnivariateDistribution}
end


function length(x::RiskParameterPriors)
  params = sum([length(x.sparks);
                length(x.susceptibility);
                length(x.transmissibility);
                length(x.infectivity)])
  if typeof(x) in [SEIR_RiskParameterPriors;
                   SEI_RiskParameterPriors]
    params += length(x.latency)
  end
  if typeof(x) in [SEIR_RiskParameterPriors;
                   SIR_RiskParameterPriors]
    params += length(x.removal)
  end
  return params
end
