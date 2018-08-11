function loglikelihood(rp::RiskParameters{T},
                       rf::RiskFunctions{T},
                       events::Events{T},
                       pop::DataFrame;
                       loglikelihood_output::Bool=true,
                       transmission_network_output::Bool=true) where T <: EpidemicModel
  # Initialize
  ll = 0.
  # update_queue = Int64[] # See TODO below
  # We can use the Simulation struct to recreate epidemic
  s = Simulation(pop, rf, rp)
  s.events = events
  # Event order
  event_order = sortperm(events[_state_progressions[T][2:end]][:])
  local last_event
  for i = 1:length(event_order)
    id, state_index = Tuple(CartesianIndices((events.individuals,
                                        length(_state_progressions[T][2:end])))[event_order[i]])
    new_state = _state_progressions[T][state_index+1]
    time = s.events[new_state][id]
    isnan(time) && break
    event = Event{T}(time, id, new_state)
    if loglikelihood_output
      # Get event rate totals
      rate_total = sum(sum(s.event_rates[new_state]) for state in _state_progressions[T][2:end])
      @logmsg LogLevel(-5000) "Event rate total $i = $(round(rate_total, 3))"
      if i > 1
        # Calculate length of inter-event period (ignore for 1st event)
        ΔT = event.time - last_event.time
        # Add event occurence contribution to loglikelihood
        ll += log(rate_total) - rate_total * ΔT
        # Stop log likelihood calculation anytime the loglikelihood goes to -Inf
        if ll == -Inf
          @debug "Event $i resulted in a -Inf loglikelihood (transition of individual $id at t = $(round(time, 3)) to state $new_state)"
          break
        end
        # Get the individual rate associated with the event that occurred
        event_rate = s.event_rates[event.new_state][event.individual]
        # Add the specific event contribution to loglikelihood
        ll += log(event_rate/rate_total)
      end
    end
    # Updates
    last_event = event
    update!(s.disease_states, event)
    if transmission_network_output
      update!(s.transmission_network,
              generate(Transmission,
                       s.transmission_rates,
                       event))
    end
    update!(s.transmission_rates,
            s.event_rates,
            event,
            s.disease_states,
            s.population,
            s.risk_functions,
            s.risk_parameters)
  end
  if transmission_network_output && loglikelihood_output
    return ll, s.transmission_network
  elseif transmission_network_output & !loglikelihood_output
    return s.transmission_network
  elseif !transmission_network_output & loglikelihood_output
    return ll
  else
    @error "You must return either a transmission_network proposal or loglikelihood"
  end
end
