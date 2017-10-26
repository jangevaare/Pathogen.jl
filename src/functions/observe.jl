"""
observe(events::SEIR_Events,
        delay_infected::UnivariateDistribution,
        delay_removed::UnivariateDistribution,
        force = false::Bool)

Make event time observations from a simulation. Force option ensures all
infections are observed.
"""
function observe(events::SEIR_Events,
                 delay_infected::UnivariateDistribution,
                 delay_removed::UnivariateDistribution,
                 force = false::Bool)
  infected = fill(NaN, events.individuals)
  removed = fill(NaN, events.individuals)
  if force
    @simd for i in find(.!isnan.(events.infected))
      infection_delay_ub = events.removed[i] - events.infected[i]
      infected[i] = events.infected[i] + rand(Truncated(delay_infected, 0., infection_delay_ub))
      removed[i] = events.removed[i] + rand(delay_removed)
    end
  else
    @simd for i in find(.!isnan.(events.infected))
      infection_delay = rand(delay_infected)
      if isnan(events.removed[i]) || infection_delay + events.infected[i] < events.removed[i]
        infected[i] = events.infected[i] + infection_delay
        removed[i] = events.removed[i] + rand(delay_removed)
      end
    end
  end
  return SEIR_EventObservations(infected, removed)
end


"""
observe(events::SIR_Events,
        delay_infected::UnivariateDistribution,
        delay_removed::UnivariateDistribution,
        force = false::Bool)

Make event time observations from a simulation. Force option ensures all
infections are observed.
"""
function observe(events::SIR_Events,
                 delay_infected::UnivariateDistribution,
                 delay_removed::UnivariateDistribution,
                 force = false::Bool)
  infected = fill(NaN, events.individuals)
  removed = fill(NaN, events.individuals)
  if force
    @simd for i in find(!isnan.(events.infected))
      infection_delay_ub = events.removed[i] - events.infected[i]
      infected[i] = events.infected[i] + rand(Truncated(delay_infected, 0., infection_delay_ub))
      removed[i] = events.removed[i] + rand(delay_removed)
    end
  else
    @simd for i in find(!isnan.(events.infected))
      infection_delay = rand(delay_infected)
      if isnan(events.removed[i]) || infection_delay + events.infected[i] < events.removed[i]
        infected[i] = events.infected[i] + infection_delay
        removed[i] = events.removed[i] + rand(delay_removed)
      end
    end
  end
  return SIR_EventObservations(infected, removed)
end


"""
observe(events::SEI_Events,
        delay_infected::UnivariateDistribution)

Make event time observations from a simulation. Force option ensures all
infections are observed.
"""
function observe(events::SEI_Events,
                 delay_infected::UnivariateDistribution)
  infected = fill(NaN, events.individuals)
  @simd for i in find(!isnan.(events.infected))
    infected[i] = events.infected[i] + rand(delay_infected)
  end
  return SEI_EventObservations(infected)
end


"""
observe(events::SI_Events,
        delay_infected::UnivariateDistribution)

Make event time observations from a simulation. Force option ensures all
infections are observed.
"""
function observe(events::SI_Events,
                 delay_infected::UnivariateDistribution)
  infected = fill(NaN, events.individuals)
  @simd for i in find(!isnan.(events.infected))
    infected[i] = events.infected[i] + rand(delay_infected)
  end
  return SI_EventObservations(infected)
end
