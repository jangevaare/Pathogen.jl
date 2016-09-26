"""
Prior distributions vectors for the `RiskParameters`
"""
type RiskParameterPriors
  sparks::Vector{UnivariateDistribution}
  susceptibility::Vector{UnivariateDistribution}
  transmissibility::Vector{UnivariateDistribution}
  infectivity::Vector{UnivariateDistribution}
  latency::Vector{UnivariateDistribution}
  removal::Vector{UnivariateDistribution}
end


function rand(riskpriors::RiskParameterPriors)
  sparks = Float64[]
  susceptibility = Float64[]
  transmissibility = Float64[]
  infectivity = Float64[]
  latency = Float64[]
  removal = Float64[]

  for i = 1:length(riskpriors.sparks)
    push!(sparks, rand(riskpriors.sparks[i]))
  end

  for i = 1:length(riskpriors.susceptibility)
    push!(susceptibility, rand(riskpriors.susceptibility[i]))
  end

  for i = 1:length(riskpriors.transmissability)
    push!(transmissability, rand(riskpriors.transmissability[i]))
  end

  for i = 1:length(riskpriors.infectivity)
    push!(infectivity, rand(riskpriors.infectivity[i]))
  end

  for i = 1:length(riskpriors.latency)
    push!(latency, rand(riskpriors.latency[i]))
  end

  for i = 1:length(riskpriors.removal)
    push!(removal, rand(riskpriors.removal[i]))
  end

  return RiskParameters(sparks,
                        susceptibility,
                        transmissability,
                        infectivity,
                        latency,
                        removal)
end


function logprior(riskparams::RiskParameters,
                  riskpriors::RiskParameterPriors)
  lp = 0.
  for i = 1:length(riskparams.sparks)
    lp += loglikelihood(riskpriors.sparks[i], riskparams.sparks[i])
  end
  for i = 1:length(riskparams.susceptibility)
    lp += loglikelihood(riskpriors.susceptibility[i], riskparams.susceptibility[i])
  end
  for i = 1:length(riskparams.transmissability)
    lp += loglikelihood(riskpriors.transmissability[i], riskparams.transmissability[i])
  end
  for i = 1:length(riskparams.infectivity)
    lp += loglikelihood(riskpriors.infectivity[i], riskparams.infectivity[i])
  end
  for i = 1:length(riskparams.latency)
    lp += loglikelihood(riskpriors.latency[i], riskparams.latency[i])
  end
  for i = 1:length(riskparams.removal)
    lp += loglikelihood(riskpriors.removal[i], riskparams.removal[i])
  end
  return lp
end


"""
Priors for event times
"""
type EventPriors
  exposed::Vector{Nullable{UnivariateDistribution}}
  infected::Vector{Nullable{UnivariateDistribution}}
  removed::Vector{Nullable{UnivariateDistribution}}
end


"""
Calculate log priors
"""
function logprior(events::Events, priors::EventPriors)
  lp = 0.
  for i = 1:length(events.exposed)
    if !isnull(priors.exposed[i])
      lp += loglikelihood(get(priors.exposed[i]), events.exposed[i])
    end
  end
  for i = 1:length(events.infected)
    if !isnull(priors.infected[i])
      lp += loglikelihood(get(priors.infected[i]), events.infected[i])
    end
  end
  for i = 1:length(events.removed)
    if !isnull(priors.removed[i])
      lp += loglikelihood(get(priors.removed[i]), events.removed[i])
    end
  end
  return lp
end


function rand(eventpriors::EventPriors)
  exposed = fill(NaN, length(eventpriors.exposed))
  infected = fill(NaN, length(eventpriors.infected))
  removed = fill(NaN, length(eventpriors.removed))
  for i = 1:length(exposed)
    if !isnull(eventpriors.exposed[i])
      exposed[i] = rand(eventpriors.exposed[i])
    end
  end
  for i = 1:length(infected)
    if !isnull(eventpriors.infected[i])
      infected[i] = rand(eventpriors.infected[i])
    end
  end
  for i = 1:length(removed)
    if !isnull(eventpriors.removed[i])
      removed[i] = rand(eventpriors.removed[i])
    end
  end
  return Events(exposed, infected, removed)
end
