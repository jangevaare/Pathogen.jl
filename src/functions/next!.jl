function next!(s::Simulation{T}) where T <: EpidemicModel
  event = generate(Event{T}, s.event_rates, s.time)
  if event.time < Inf
    update!(s.events, event)
    update!(s.disease_states, event)
    # Checks if a _new_transmission within function
    # Could make this into a single function for conveinence (also appears in Simulation.jl)
    update!(s.transmission_network,
            generate(Transmission,
                     s.transmission_rates,
                     event))
    # Update `EventRates` before `TransmissionRates`
    update!(s.event_rates,
            s.transmission_rates,
            event,
            s.disease_states,
            s.population,
            s.risk_functions,
            s.risk_parameters)
    update!(s.transmission_rates,
            event,
            s.disease_states,
            s.population,
            s.risk_functions,
            s.risk_parameters)
    # Count the iteration as having occurred.
    s.iterations += 1
  end
  # Update simulation time
  s.time = event.time
  # Return updated Simulation object
  return s
end

function _accept(lp1::Float64,
                 lp2::Float64)
  return rand() <= exp(lp1 - lp2)
end

function next!(mc::MarkovChain{T},
               mcmc::MCMC{T},
               Σ::Array{Float64, 2},
               σ::Float64) where T <: EpidemicModel
  # Initialize
  new_events = mc.events[end]
  new_events_array = new_events[_state_progressions[T][2:end]]
  new_params = mc.risk_parameters[end]
  new_network = mc.network[end]
  new_lposterior = mc.log_posterior[end]
  # Randomize event time augmentation
  event_indices = find(.!isnan.(new_events_array[:]))
  augmentation_order = sample(event_indices, length(event_indices), replace=false)
  for i = 1:(length(event_indices) + 1)
    if i <= length(event_indices)
      id, state_index = ind2sub((new_events.individuals, length(_state_progressions[T][2:end])),
                                  augmentation_order[i])
      new_state = _state_progressions[T][state_index+1]
      time = new_events[new_state][id]
      proposed_event = generate(Event{T},
                                Event{T}(time, id, new_state),
                                σ,
                                mcmc.event_extents,
                                mcmc.event_observations,
                                new_events)
      proposed_events_array = reshape([new_events_array[1:(augmentation_order[i]-1)]
                                       proposed_event.time
                                       new_events_array[(augmentation_order[i]+1):end]],
                                      size(new_events_array))
      proposed_events = Events{T}(proposed_events_array)
      proposed_params = new_params
    else
      # Propose new risk parameters
      proposed_events = new_events
      proposed_events_array = new_events_array
      proposed_params = generate(RiskParameters{T}, new_params, Σ)
    end
    proposed_lprior = logpriors(proposed_params, mcmc.risk_priors)
    if proposed_lprior > -Inf
      proposed_lliklihood, proposed_network = loglikelihood(proposed_params,
                                                            mcmc.risk_functions,
                                                            proposed_events,
                                                            mcmc.population)
      if proposed_lliklihood > -Inf
        proposed_lposterior = proposed_lprior + proposed_lliklihood
      else
        proposed_lposterior = -Inf
      end
    else
      proposed_lposterior = -Inf
    end
    if _accept(proposed_lposterior, new_lposterior)
      new_events = proposed_events
      new_events_array = proposed_events_array
      new_params = proposed_params
      new_network = proposed_network
      new_lposterior = proposed_lposterior
    end
  end
  mc.iterations += 1
  push!(mc.events, new_events)
  push!(mc.network, new_network)
  push!(mc.risk_parameters, new_params)
  push!(mc.log_posterior, new_lposterior)
  return mc
end

function next!(mcmc::MCMC{T},
               Σ::Array{Float64, 2},
               σ::Float64) where T <: EpidemicModel
  @simd for mc in mcmc.markov_chains
    next!(mc, mcmc, Σ, σ)
  end
  return mcmc
end
