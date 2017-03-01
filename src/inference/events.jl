"""
generate_events(observations::EventObservations,
                exposureextent::Float64,
                infectionextent::Float64,
                removalextent::Float64)

Generate some initial `Events` from `EventObservations`
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
      removed_ub = observations.removed[i]
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
propose(i::Int64,
        j::Int64,
        events::Events,
        network::Network,
        observations::EventObservations,
        variance::Float64)

Event time augmentation for a single event time
"""
function propose(i::Int64,
                 j::Int64,
                 events::Events,
                 network::Network,
                 observations::EventObservations,
                 variance::Float64)
  exposed = copy(events.exposed)
  infected = copy(events.infected)
  removed = copy(events.removed)
  pathfrom = pathwayfrom(i, network, 2)
  pathto = pathwayto(i, network, 2)
  # Exposure time
  if j == 1
    if length(pathto) > 1
      exposure_lb = maximum([infected[pathto[2]]; -Inf])
      exposure_ub = minimum([infected[i]; removed[pathto[2]]])
      exposed[i] = rand(TruncatedNormal(exposed[i], variance, exposure_lb, exposure_ub))
    else
      exposure_lb = -Inf
      exposure_ub = infected[i]
      exposed[i] = rand(TruncatedNormal(exposed[i], variance, exposure_lb, exposure_ub))
    end
  # Infection time
  elseif j == 2
    infection_lb = exposed[i]
    infection_ub = minimum([exposed[pathfrom[2:end]]; observations.infected[i]; Inf])
    infected[i] = rand(TruncatedNormal(infected[i], variance, infection_lb, infection_ub))
  # Removal time
  elseif j == 3
    removal_lb = maximum([exposed[pathfrom[2:end]]; observations.infected[i]; -Inf])
    removal_ub = Inf
    removed[i] = rand(TruncatedNormal(removed[i], variance, removal_lb, removal_ub))
  end
  return Events(exposed,
                infected,
                removed)
end
