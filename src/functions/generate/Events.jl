function generate(::Type{Events},
                  rp::RiskParameters{T},
                  rf::RiskFunctions{T},
                  pop::Population,
                  obs::EventObservations{T},
                  extents::EventExtents{T}) where T <: EpidemicModel
  events = Events{T}(obs.individuals)
  exposed_state = T in [SEIR; SEI]
  removed_state = T in [SEIR; SIR]
  for i = 1:obs.individuals
    if obs.infection[i] == -Inf
      update!(events, Event{T}(-Inf, i, State_I))
      if exposed_state
        update!(events, Event{T}(-Inf, i, State_E))
      end
      if removed_state && obs.removal[i] == -Inf
        update!(events, Event{T}(-Inf, i, State_R))
      end
    elseif !isnan(obs.infection[i])
      i_lb = obs.infection[i] - extents.infection
      i_ub = obs.infection[i]
      # Generate infection time as a uniformly distributed random variable
      # Event time is bounded by specified infection extent, and the actual infection observation
      i_time = rand(Uniform(i_lb, i_ub))
      update!(events, Event{T}(i_time, i, State_I))
      if exposed_state
        # Calculate individual specific latency rate
        λ = rf.latency(rp.latency, pop, i)
        # Generate the exposure time as a (truncated) exponentially distributed random variable
        # The latent period can not exceed the specified latent period extent
        e_time = i_time - rand(Truncated(Exponential(1.0/λ), 0.0, extents.exposure))
        # Update the `Events` record
        update!(events, Event{T}(e_time, i, State_E))
      end
      if removed_state && !isnan(obs.removal[i])
        # Calculate individual specific removal rate
        λ = rf.removal(rp.removal, pop, i)
        # Find the lower bound for removal time...
        # Removal event must be AFTER...
        # - The infection observation
        # - (Removal observation time - max removal observation delay)
        r_lb = maximum([obs.infection[i]; obs.removal[i] - extents.removal])
        # Find the upper bound for removal time...
        # Removal event must be prior to removal observation
        r_ub = obs.removal[i]
        # Generate the removal time as a (truncated) exponetially distributed random variable
        r_time = rand(Truncated(Exponential(1.0/λ), r_lb, r_ub))
        # Update the `Events` record
        update!(events, Event{T}(r_time, i, State_R))
      end
    end
  end
  return events
end

function generate(::Type{Events},
                  rp::RiskParameters{T},
                  mcmc::MCMC{T}) where T <: EpidemicModel
  return generate(Events,
                  rp,
                  mcmc.risk_functions,
                  mcmc.population,
                  mcmc.event_observations,
                  mcmc.event_extents)
end
