function generate(::Type{Events},
                  obs::EventObservations{S},
                  extents::EventExtents{S}) where S <: DiseaseStateSequence
  events = Events{S}(individuals(obs))
  exposed_state = S in [SEIR; SEI]
  removed_state = S in [SEIR; SIR]
  for i = 1:individuals(obs)
    if obs.infection[i] == -Inf
      update!(events, Event{S}(-Inf, i, State_I))
      if exposed_state
        update!(events, Event{S}(-Inf, i, State_E))
      end
      if removed_state && obs.removal[i] == -Inf
        update!(events, Event{S}(-Inf, i, State_R))
      end
    elseif !isnan(obs.infection[i])
      i_lb = maximum([obs.infection[i] - extents.infection[2]])
      i_ub = obs.infection[i] - extents.infection[1]
      i_time = rand(Uniform(i_lb, i_ub))
      update!(events, Event{S}(i_time, i, State_I))
      if exposed_state
        e_lb = maximum([i_time - extents.exposure[2]])
        e_ub = i_time - extents.exposure[1]
        e_time = rand(Uniform(e_lb, e_ub))
        update!(events, Event{S}(e_time, i, State_E))
      end
      if removed_state && !isnan(obs.removal[i])
        r_lb = maximum([obs.infection[i]; obs.removal[i] - extents.removal[2]])
        r_ub = obs.removal[i] - extents.removal[1]
        r_time = rand(Uniform(r_lb, r_ub))
        update!(events, Event{S}(r_time, i, State_R))
      end
    end
  end
  return events
end

function generate(::Type{Events}, mcmc::MCMC)
  return generate(Events, mcmc.event_observations, mcmc.event_extents)
end