function observe(events::Events{T},
                 delay_infection::UnivariateDistribution,
                 delay_removal::UnivariateDistribution;
                 force = false::Bool) where T <: Union{SEIR, SIR}
  infection = fill(NaN, events.individuals)
  removal = fill(NaN, events.individuals)
  if force
    @simd for i in find(.!isnan.(events.infection))
      if isnan(events.removal[i])
        infection[i] = events.infection[i] + rand(delay_infection)
      else
        infection_delay_ub = events.removal[i] - events.infection[i]
        infection[i] = events.infection[i] + rand(Truncated(delay_infection, 0., infection_delay_ub))
        removal[i] = events.removal[i] + rand(delay_removal)
      end
    end
  else
    @simd for i in find(.!isnan.(events.infection))
      infection_delay = rand(delay_infection)
      if isnan(events.removal[i]) || infection_delay + events.infection[i] < events.removal[i]
        infection[i] = events.infection[i] + infection_delay
        if !isnan(events.removal[i])
          removal[i] = events.removal[i] + rand(delay_removal)
        end
      end
    end
  end
  return EventObservations{T}(infection, removal)
end

function observe(sim::Simulation{T},
                 delay_infection::UnivariateDistribution,
                 delay_removal::UnivariateDistribution;
                 force = false::Bool) where T <: Union{SEIR, SIR}
  return observe(sim.events, delay_infection, delay_removal, force=force)
end

function observe(events::Events{T},
                 delay_infection::UnivariateDistribution) where T <: Union{SEI, SI}
  infection = fill(NaN, events.individuals)
  @simd for i in find(.!isnan.(events.infection))
    infection[i] = events.infection[i] + rand(delay_infection)
  end
  return EventObservations{T}(infection)
end

function observe(sim::Simulation{T},
                 delay_infection::UnivariateDistribution) where T <: Union{SEI, SI}
  return observe(sim.events, delay_infection)
end
