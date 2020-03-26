function generate(::Type{Events},
                  obs::EventObservations{S},
                  extents::EventExtents{S}) where S <: DiseaseStateSequence
  events = Events{S}(obs.individuals)
  exposed_state = S in [SEIR; SEI]
  removed_state = S in [SEIR; SIR]
  for i = 1:obs.individuals
    if obs.infection[i] == -Inf
      update!(events, Event{S}(-Inf, i, State_I))
      if exposed_state
        update!(events, Event{S}(-Inf, i, State_E))
      end
      if removed_state && obs.removal[i] == -Inf
        update!(events, Event{S}(-Inf, i, State_R))
      end
    elseif !isnan(obs.infection[i])
      i_lb = obs.infection[i] - extents.infection
      i_ub = obs.infection[i]
      i_time = rand(Uniform(i_lb, i_ub))
      update!(events, Event{S}(i_time, i, State_I))
      if exposed_state
        e_lb = i_time - extents.exposure
        e_ub = i_time
        e_time = rand(Uniform(e_lb, e_ub))
        update!(events, Event{S}(e_time, i, State_E))
      end
      if removed_state && !isnan(obs.removal[i])
        r_lb = maximum([obs.infection[i]; obs.removal[i] - extents.removal])
        r_ub = obs.removal[i]
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