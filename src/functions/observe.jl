function observe(events::Vector{Event{T}},
                 delay_infection::UnivariateDistribution,
                 delay_removal::UnivariateDistribution,
                 n_ids::Int64,
                 force = false::Bool) where T <: Union{SEIR, SIR}
  infection_obs = fill(NaN, n_ids)
  removal_obs = fill(NaN, n_ids)
  removal = fill(NaN, n_ids)
  event_order = sortperm(_time.(events))
  for i in reverse(event_order)
    if events[i].new_state == State_R
      removal[events[i].individual] = events[i].time
      time_obs = events[i].time + rand(delay_removal)
      removal_obs[events[i].individual] = time_obs
    elseif events[i].new_state == State_I
      if !isnan(removal[events[i].individual])
        if force
          ub = removal[events[i].individual] - events[i].time
          time_obs = events[i].time + rand(Truncated(delay_infection, 0, ub))
          infection_obs[events[i].individual] = time_obs
        else
          time_obs = events[i].time + rand(delay_infection)
          if time_obs < removal[events[i].individual]
            infection_obs[events[i].individual] = time_obs
          else
            infection_obs[events[i].individual] = NaN
            removal_obs[events[i].individual] = NaN
          end
        end
      end
    end
  end
  return EventObservations{T}(infection_obs, removal_obs)
end

function observe(events::Vector{Event{T}},
                 delay_infection::UnivariateDistribution,
                 n_ids::Int64) where T <: Union{SEI, SI}
  infection_obs = fill(NaN, n_ids)
  event_order = sortperm(_time.(events))
  for i in event_order
    if events[i].new_state == State_I
      time_obs = events[i].time + rand(delay_infection)
      infection_obs[events[i].individual] = time_obs
    end
  end
  return EventObservations{T}(infection_obs)
end

function observe(sim::Simulation{T},
                 delay_infection::UnivariateDistribution) where T <: Union{SEI, SI}
  return observe(sim.events, delay_infection, size(sim.population, 1))
end

function observe(sim::Simulation{T},
                 delay_infection::UnivariateDistribution,
                 delay_removal::UnivariateDistribution,
                 force = false::Bool) where T <: Union{SEIR, SIR}
  return observe(sim.events, delay_infection, delay_removal, size(sim.population, 1), force)
end
