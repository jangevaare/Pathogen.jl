function loglikelihood(rp::RiskParameters{T},
                       rf::RiskFunctions{T},
                       events::Events{T},
                       pop::DataFrame;
                       debug_level::Int64=0) where T <: EpidemicModel
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
    individual, state_index = ind2sub((events.individuals, length(_state_progressions[T][2:end])),
                                      event_order[i])
    new_state = _state_progressions[T][state_index+1]
    time = s.events[new_state][individual]
    isnan(time) && break
    event = Event{T}(time, individual, new_state)
    # Get event rate totals
    rate_total = sum(sum(s.event_rates[new_state]) for state in _state_progressions[T][2:end])
    if debug_level >= 4
      println("loglikelihood: event rate total $i = $(round(rate_total, 3))")
    end
    if i > 1
      # Calculate length of inter-event period (ignore for 1st event)
      ΔT = event.time - last_event.time
      # Add event occurence contribution to loglikelihood
      ll += log(rate_total) - rate_total * ΔT
      # Stop log likelihood calculation anytime the loglikelihood goes to -Inf
      if ll == -Inf
        if debug_level >= 2
          println("loglikelihood: event $i resulted in a -Inf loglikelihood (transition of individual $individual at t = $(round(time, 3)) to state $new_state)")
        end
        break
      end
      # Get the individual rate associated with the event that occurred
      event_rate = s.event_rates[event.new_state][event.individual]
      # Add the specific event contribution to loglikelihood
      ll += log(event_rate/rate_total)
    end

    # Updates
    last_event = event
    update!(s.disease_states, event)
    update!(s.transmission_network,
            generate(Transmission,
                     s.transmission_rates,
                     event,
                     debug_level = debug_level))
    update!(s.transmission_rates,
            s.event_rates,
            event,
            s.disease_states,
            s.population,
            s.risk_functions,
            s.risk_parameters,
            debug_level = debug_level)

    # TODO Consider using a queue to eliminate sensitivity to simulataneous events
    # Add to queue
    # push!(update_queue, i)
    # if i > 1 && ΔT > 0.0
    #   for j in update_queue
    #     queued_event = events[event_order[j]]
    #     update!(s.events, queued_event)
    #     update!(s.disease_states, queued_event)
    #     update!(s.transmission_network,
    #             generate(Transmission,
    #                      s.transmission_rates,
    #                      queued_event))
    #     # Will need to break these into two parts to make work
    #     update!(s.event_rates,
    #             s.transmission_rates,
    #             queued_event,
    #             s.disease_states,
    #             s.population,
    #             s.risk_functions,
    #             s.risk_parameters)
    #     update!(s.transmission_rates,
    #             queued_event,
    #             s.disease_states,
    #             s.population,
    #             s.risk_functions,
    #             s.risk_parameters)
    #   end
    #   update_queue = Int64[]
    # end
  end
  return ll, s.transmission_network
end
