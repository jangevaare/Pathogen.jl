"""
generate_events(observations::SEIR_EventObservations,
                exposureextent::Float64,
                infectionextent::Float64,
                removalextent::Float64)

Generate some initial `SEIR_Events` from `SEIR_EventObservations`
"""
function generate_events(observations::SEIR_EventObservations,
                         exposureextent::Float64,
                         infectionextent::Float64,
                         removalextent::Float64)
  removed = fill(NaN, observations.individuals)
  infected = fill(NaN, observations.individuals)
  exposed = fill(NaN, observations.individuals)
  for i = 1:observations.individuals
    if !isnan(observations.removed[i])
      removed_lb = maximum([observations.infected[i];
                            observations.removed[i]-removalextent])
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
  return SEIR_Events(exposed,
                     infected,
                     removed)
end


"""
generate_events(observations::SIR_EventObservations,
                exposureextent::Float64,
                infectionextent::Float64,
                removalextent::Float64)

Generate some initial `SIR_Events` from `SIR_EventObservations`
"""
function generate_events(observations::SIR_EventObservations,
                         infectionextent::Float64,
                         removalextent::Float64)
  removed = fill(NaN, observations.individuals)
  infected = fill(NaN, observations.individuals)
  for i = 1:observations.individuals
    if !isnan(observations.removed[i])
      removed_lb = maximum([observations.infected[i];
                            observations.removed[i]-removalextent])
      removed_ub = observations.removed[i]
      removed[i] = rand(Uniform(removed_lb, removed_ub))
    end
    if !isnan(observations.infected[i])
      infected_lb = observations.infected[i] - infectionextent
      infected_ub = observations.infected[i]
      infected[i] = rand(Uniform(infected_lb, infected_ub))
    end
  end
  return SIR_Events(infected,
                    removed)
end


"""
generate_events(observations::SEI_EventObservations,
                exposureextent::Float64,
                infectionextent::Float64,
                removalextent::Float64)

Generate some initial `SEI_Events` from `SEI_EventObservations`
"""
function generate_events(observations::SEI_EventObservations,
                         exposureextent::Float64,
                         infectionextent::Float64)
  infected = fill(NaN, observations.individuals)
  exposed = fill(NaN, observations.individuals)
  for i = 1:observations.individuals
    if !isnan(observations.infected[i])
      infected_lb = observations.infected[i] - infectionextent
      infected_ub = observations.infected[i]
      infected[i] = rand(Uniform(infected_lb, infected_ub))
      exposed_lb = infected[i] - exposureextent
      exposed_ub = infected[i]
      exposed[i] = rand(Uniform(exposed_lb, exposed_ub))
    end
  end
  return SEI_Events(exposed,
                    infected)
end


"""
generate_events(observations::SI_EventObservations,
                exposureextent::Float64,
                infectionextent::Float64,
                removalextent::Float64)

Generate some initial `SI_Events` from `SI_EventObservations`
"""
function generate_events(observations::SIR_EventObservations,
                         infectionextent::Float64)
  infected = fill(NaN, observations.individuals)
  for i = 1:observations.individuals
    if !isnan(observations.infected[i])
      infected_lb = observations.infected[i] - infectionextent
      infected_ub = observations.infected[i]
      infected[i] = rand(Uniform(infected_lb, infected_ub))
    end
  end
  return SIR_Events(infected)
end
