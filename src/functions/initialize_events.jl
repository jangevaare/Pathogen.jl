"""
initialize_events(observations::SEIR_EventObservations,
                  event_extents::SEIR_EventExtents)

Generate some initial `SEIR_Events` from `SEIR_EventObservations`
"""
function initialize_events(observations::SEIR_EventObservations,
                           event_extents::SEIR_EventExtents)
  removed = fill(NaN, observations.individuals)
  infected = fill(NaN, observations.individuals)
  exposed = fill(NaN, observations.individuals)
  for i = 1:observations.individuals
    if !isnan(observations.removed[i])
      removed_lb = maximum([observations.infected[i];
                            observations.removed[i]-event_extents.removal])
      removed_ub = observations.removed[i]
      removed[i] = rand(Uniform(removed_lb, removed_ub))
    end
    if !isnan(observations.infected[i])
      infected_lb = observations.infected[i] - event_extents.infection
      infected_ub = observations.infected[i]
      infected[i] = rand(Uniform(infected_lb, infected_ub))
      exposed_lb = infected[i] - event_extents.exposure
      exposed_ub = infected[i]
      exposed[i] = rand(Uniform(exposed_lb, exposed_ub))
    end
  end
  return SEIR_Events(exposed,
                     infected,
                     removed)
end


"""
initialize_events(observations::SIR_EventObservations,
                  event_extents::SIR_EventExtents)

Generate some initial `SIR_Events` from `SIR_EventObservations`
"""
function initialize_events(observations::SIR_EventObservations,
                           event_extents::SIR_EventExtents)
  removed = fill(NaN, observations.individuals)
  infected = fill(NaN, observations.individuals)
  for i = 1:observations.individuals
    if !isnan(observations.removed[i])
      removed_lb = maximum([observations.infected[i];
                            observations.removed[i]-event_extents.removal])
      removed_ub = observations.removed[i]
      removed[i] = rand(Uniform(removed_lb, removed_ub))
    end
    if !isnan(observations.infected[i])
      infected_lb = observations.infected[i] - event_extents.infection
      infected_ub = observations.infected[i]
      infected[i] = rand(Uniform(infected_lb, infected_ub))
    end
  end
  return SIR_Events(infected,
                    removed)
end


"""
initialize_events(observations::SEI_EventObservations,
                  event_extents::SEI_EventExtents)

Generate some initial `SEI_Events` from `SEI_EventObservations`
"""
function initialize_events(observations::SEI_EventObservations,
                         event_extents::SEI_EventExtents)
  infected = fill(NaN, observations.individuals)
  exposed = fill(NaN, observations.individuals)
  for i = 1:observations.individuals
    if !isnan(observations.infected[i])
      infected_lb = observations.infected[i] - event_extents.infection
      infected_ub = observations.infected[i]
      infected[i] = rand(Uniform(infected_lb, infected_ub))
      exposed_lb = infected[i] - event_extents.exposure
      exposed_ub = infected[i]
      exposed[i] = rand(Uniform(exposed_lb, exposed_ub))
    end
  end
  return SEI_Events(exposed,
                    infected)
end


"""
initialize_events(observations::SI_EventObservations,
                  event_extents::SI_EventExtents)

Generate some initial `SI_Events` from `SI_EventObservations`
"""
function initialize_events(observations::SI_EventObservations,
                           event_extents::SI_EventExtents)
  infected = fill(NaN, observations.individuals)
  for i = 1:observations.individuals
    if !isnan(observations.infected[i])
      infected_lb = observations.infected[i] - event_extents.infection
      infected_ub = observations.infected[i]
      infected[i] = rand(Uniform(infected_lb, infected_ub))
    end
  end
  return SI_Events(infected)
end
