function observe(events::Events{T},
                 delay_infection::UnivariateDistribution,
                 delay_removal::UnivariateDistribution;
                 force::Bool = false) where T <: Union{SEIR, SIR}
  infection = fill(NaN, events.individuals)
  removal = fill(NaN, events.individuals)
  if force
    for i in findall(.!isnan.(events.infection))
      if events.infection[i] == -Inf
        infection[i] = -Inf
        if events.removal[i] == -Inf
          removal[i] = -Inf
        end
      elseif isnan(events.removal[i])
        infection[i] = events.infection[i] + rand(delay_infection)
      elseif events.removal[i] > -Inf
        infection_delay_ub = events.removal[i] - events.infection[i]
        infection[i] = events.infection[i] +
                       rand(Truncated(delay_infection, 0.0, infection_delay_ub))
        removal[i] = events.removal[i] + rand(delay_removal)
      end
      @debug "Infection observation of i = $i at t = $(round(infection[i], digits=3)) (actual infection onset at t = $(round(events.infection[i], digits=3)))"
      @debug "Removal observation of i = $i at t = $(round(removal[i], digits=3)) (actual removal at t = $(round(events.removal[i], digits=3)))"
    end
  else
    for i in findall(.!isnan.(events.infection))
      if events.infection[i] == -Inf
        infection[i] = -Inf
        if events.removal == -Inf
          removal[i] = -Inf
        end
      else
        infection_delay = rand(delay_infection)
        if isnan(events.removal[i])
          infection[i] = events.infection[i] + infection_delay
        else
          if infection_delay + events.infection[i] < events.removal[i]
            infection[i] = events.infection[i] + infection_delay
            removal[i] = events.removal[i] + rand(delay_removal)
          end
        end
      end
      @debug "Infection observation of i = $i t = $(round(infection[i], digits=3)) (actual infection onset at t = $(round(events.infection[i], digits=3)))"
      @debug "Removal observation of i = $i at t = $(round(removal[i], digits=3)) (actual removal at t = $(round(events.removal[i], digits=3)))"
    end
  end
  return EventObservations{T}(infection, removal)
end

function observe(sim::Simulation{T},
                 delay_infection::UnivariateDistribution,
                 delay_removal::UnivariateDistribution;
                 force::Bool = false) where T <: Union{SEIR, SIR}
  return observe(sim.events, delay_infection, delay_removal, force=force)
end

function observe(events::Events{T},
                 delay_infection::UnivariateDistribution) where T <: Union{SEI, SI}
  infection = fill(NaN, events.individuals)
  @simd for i in findall(.!isnan.(events.infection))
    if events.infection[i] == -Inf
      infection[i] = -Inf
    else
      infection[i] = events.infection[i] + rand(delay_infection)
    end
    @debug "Infection observation of i = $i at t = $(round(infection[i], digits=3))) (actual infection onset at t = $(round(events.infection[i], digits=3)))"
  end
  return EventObservations{T}(infection)
end

function observe(sim::Simulation{T},
                 delay_infection::UnivariateDistribution) where T <: Union{SEI, SI}
  return observe(sim.events, delay_infection)
end
