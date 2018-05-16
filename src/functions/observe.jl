function observe(events::Events{T},
                 delay_infected::UnivariateDistribution,
                 delay_removed::UnivariateDistribution,
                 force = false::Bool) where T <: Union{SEIR, SIR}
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
  return EventObservations{T}(infected, removed)
end

function observe(events::Events{T},
                 delay_infected::UnivariateDistribution,
                 force = false::Bool) where T <: Union{SEI, SI}
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
  return EventObservations{T}(infected)
end
