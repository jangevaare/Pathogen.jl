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
function generate_events(observations::EventObservations,
                         exposureextent::Float64,
                         infectionextent::Float64,
                         removalextent::Float64)
  removed = fill(NaN, observations.individuals)
  infected = fill(NaN, observations.individuals)
  exposed = fill(NaN, observations.individuals)
  for i = 1:observations.individuals
    if !isnan(observations.removed[i])
      removed_lb = maximum([observations.infected[i]; observations.removed[i]-removalextent])
      removed_ub = observations.removed[i])
      removed[i] = rand(Uniform(removed_lb, removed_ub))
    end
    if !isnan(observations.infected[i])
      infected_lb = observations.infected[i] - infectionextent
      infected_ub = observations.infected[i]
      infected[i] = rand(Uniform(infected_lb, infected_ub))
      exposed_lb = infected[i] - exposureextent
      exposed_ub = infected[i]
      exposed[i] = rand(Uniform(exposed_lb, exposed_ub))
    end
  end
  return Events(exposed,
                infected,
                removed)
end


"""
Calculate log priors
"""
function logprior(priors::EventPriors,
                  events::Events)
  lp = 0.
  @simd for i = 1:length(events.exposed)
    if !isnull(priors.exposed[i])
      lp += loglikelihood(get(priors.exposed[i]), [events.exposed[i]])
    end
  end
  @simd for i = 1:length(events.infected)
    if !isnull(priors.infected[i])
      lp += loglikelihood(get(priors.infected[i]), [events.infected[i]])
    end
  end
  @simd for i = 1:length(events.removed)
    if !isnull(priors.removed[i])
      lp += loglikelihood(get(priors.removed[i]), [events.removed[i]])
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
  exposed = events.exposed
  infected = events.infected
  removed = events.removed
  for i in individuals
    pathfrom = pathwayfrom(i, network, 2)
    pathto = pathwayto(i, network, 2)
    # Exposure time
    if !isnan(exposed[i])
      if length(pathto) > 1
        exposure_lb = infected[pathto[2]]
        if isnan(removed[pathto[2]])
          exposure_ub = infected[i]
        else
          exposure_ub = minimum([infected[i]; exposed[pathfrom[2:end]]; removed[pathto[2]]])
        end
      else
        exposure_lb = 0.
        exposure_ub = infected[i]
      end
      exposed[i] = rand(Truncated(get(eventpriors.exposed[i]), exposure_lb, exposure_ub))
    end
    # Infection time
    if !isnan(infected[i])
      infection_lb = exposed[i]
      infection_ub = minimum([exposed[pathfrom[2:end]]; Inf])
      infected[i] = rand(Truncated(get(eventpriors.infected[i]), infection_lb, infection_ub))
    end
    # Removal time
    if !isnan(removed[i])
      removal_lb = maximum([exposed[pathfrom[2:end]]; infected[i]])
      removal_ub = Inf
      removed[i] = rand(Truncated(get(eventpriors.removed[i]), removal_lb, removal_ub))
    end
  end
  return Events(exposed, infected, removed)
end


"""
Event time augmentation for a single event time
"""
function propose(events::Events,
                 network::Network,
                 variance::Float64)
  individuals = events.individuals
  exposed = events.exposed
  infected = events.infected
  removed = events.removed
  # Randomly select an event time
  i, j = ind2sub((individuals, 3), sample(find([!isnan(exposed) !isnan(infected) !isnan(removed)]), 1)[1])
  pathfrom = pathwayfrom(i, network, 2)
  pathto = pathwayto(i, network, 2)
  # Exposure time
  if j == 1
    if length(pathto) > 1
      exposure_lb = infected[pathto[2]]
      exposure_ub = minimum([infected[i]; removed[pathto[2]]; Inf])
    else
      exposure_lb = -Inf
      exposure_ub = minimum([infected[i]; Inf])
    end
    exposed[i] = rand(TruncatedNormal(exposed[i], variance, exposure_lb, exposure_ub))
  # Infection time
  elseif j == 2
    infection_lb = exposed[i]
    infection_ub = minimum([exposed[pathfrom[2:end]]; Inf])
    infected[i] = rand(TruncatedNormal(infected[i], variance, infection_lb, infection_ub))
  # Removal time
  elseif j == 3
    removal_lb = maximum([exposed[pathfrom[2:end]]; infected[i]])
    removal_ub = Inf
    removed[i] = rand(TruncatedNormal(removed[i], variance, removal_lb, removal_ub))
  end
  return Events(exposed, infected, removed)
end


"""
Event time augmentation
"""
function propose(individuals::Vector{Int64},
                 events::Events,
                 eventpriors::EventPriors)
  exposed = events.exposed
  infected = events.infected
  removed = events.removed
  for i in individuals
    if !isnull(eventpriors.exposed[i])
      exposed[i] = rand(get(eventpriors.exposed[i]))
    end
    if !isnull(eventpriors.infected[i])
      infected[i] = rand(get(eventpriors.infected[i]))
    end
    if !isnull(eventpriors.removed[i])
      removed[i] = rand(get(eventpriors.removed[i]))
    end
  end
  return Events(exposed, infected, removed)
end
