function observe(events::Events{T},
                 delay_infected::UnivariateDistribution,
                 delay_removed::UnivariateDistribution,
                 force = false::Bool) where T <: Union{SEIR, SIR}
  infected = fill(NaN, events.individuals)
  removed = fill(NaN, events.individuals)
  if force
    for i in find(.!isnan.(events.infected))
      infection_delay_ub = events.removed[i] - events.infected[i]
      infected[i] = events.infected[i] + rand(Truncated(delay_infected, 0., infection_delay_ub))
      removed[i] = events.removed[i] + rand(delay_removed)
    end
  else
    for i in find(.!isnan.(events.infected))
      infection_delay = rand(delay_infected)
      if isnan(events.removed[i]) || infection_delay + events.infected[i] < events.removed[i]
        infected[i] = events.infected[i] + infection_delay
        removed[i] = events.removed[i] + rand(delay_removed)
      end
    end
  end
  return EventObservations{T}(infected, removed)
end

function observe(events::Events{T},
                 delay_infected::UnivariateDistribution) where T <: Union{SEI, SI}
  infected = fill(NaN, events.individuals)
  for i in find(.!isnan.(events.infected))
    infected[i] = events.infected[i] + rand(delay_infected)
  end
  return EventObservations{T}(infected)
end

function observe(sim::Simulation{T},
                 delay_infected::UnivariateDistribution) where T <: Union{SEI, SI}
  return observe(sim.events, delay_infected)
end

function observe(sim::Simulation{T},
                 delay_infected::UnivariateDistribution,
                 delay_removed::UnivariateDistribution,
                 force = false::Bool) where T <: Union{SEIR, SIR}
  return observe(sim.events, delay_infected, delay_removed, force)
end
