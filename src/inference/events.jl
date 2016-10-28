"""
Priors for event times
"""
type EventPriors
  exposed::Vector{Nullable{UnivariateDistribution}}
  infected::Vector{Nullable{UnivariateDistribution}}
  removed::Vector{Nullable{UnivariateDistribution}}
end


"""
Generate uniform `EventPriors` from `EventObservations`
"""
function generate_eventpriors(observations::EventObservations,
                              exposureextent::Real,
                              infectionextent::Real,
                              removalextent::Real)
  removed = fill(Nullable{UnivariateDistribution}(), observations.individuals)
  infected = fill(Nullable{UnivariateDistribution}(), observations.individuals)
  exposed = fill(Nullable{UnivariateDistribution}(), observations.individuals)
  for i = 1:observations.individuals
    if !isnan(observations.removed[i])
      removed_lb = maximum([observations.infected[i];
                            observations.removed[i]-removalextent])
      removed[i] = Uniform(removed_lb,
                           observations.removed[i])
    end
    if !isnan(observations.infected[i])
      infected[i] = Uniform(observations.infected[i] - infectionextent,
                            observations.infected[i])
      exposed[i] = Uniform(observations.infected[i] - infectionextent - exposureextent,
                           observations.infected[i] - infectionextent)
    end
  end
  return EventPriors(exposed,
                     infected,
                     removed)
end


"""
Calculate log priors
"""
function logprior(events::Events,
                  priors::EventPriors)
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
      exposed[i] = rand(get(eventpriors.exposed[i]))
    end
  end
  for i = 1:length(infected)
    if !isnull(eventpriors.infected[i])
      infected[i] = rand(get(eventpriors.infected[i]))
    end
  end
  for i = 1:length(removed)
    if !isnull(eventpriors.removed[i])
      removed[i] = rand(get(eventpriors.removed[i]))
    end
  end
  return Events(exposed, infected, removed)
end


"""
Event time augmentation
"""
function propose(individuals::Vector{Int64},
                 events::Events,
                 eventpriors::EventPriors,
                 network::Network)
  proposal = events
  for i in individuals
    pathfrom = pathwayfrom(i, network)
    pathto = pathwayto(i, network)
    # Exposure time
    if !isnan(events.exposed[i])
      if length(pathto) > 2
        exposure_lb = events.infected[pathto[2]]
        if isnan(events.removed[pathto[2]])
          exposure_ub = events.infected[i]
        else
          exposure_ub = minimum([events.infected[i]; events.removed[pathto[2]]])
        end
      else
        exposure_lb = 0.
        exposure_ub = events.infected[i]
      end
      proposal.exposed[i] = rand(Truncated(eventpriors.exposed[i], exposure_lb, exposure_ub))
    end
    # Infection time
    if !isnan(events.infected[i])
      infection_lb = proposal.exposed[i]
      infection_ub = minimum(proposal.exposed[pathfrom[2:end]], Inf)
      proposal.infected[i] = rand(Truncated(eventpriors.infected[i], infection_lb, infection_ub))
    end
    # Removal time
    if !isnan(events.removed[i])
      removal_lb = maximum(proposal.exposed[pathfrom[2:end]])
      removal_ub = Inf
      proposal.removed[i] = rand(Truncated(eventpriors.removed[i], removal_lb, removal_ub))
    end
  end
  return proposal
end
