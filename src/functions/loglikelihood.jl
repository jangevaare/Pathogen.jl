function loglikelihood(rp::RiskParameters{T},
                       rf::RiskFunctions{T},
                       events::Events{T},
                       pop::Population,
                       starting_states::Vector{DiseaseState};
                       loglikelihood_output::Bool=true,
                       transmission_network_output::Bool=true,
                       early_decision_value::Float64=-Inf) where T <: EpidemicModel
  # Initialize
  ll = 0.0
  # We can use the Simulation struct to recreate epidemic
  s = Simulation(pop, starting_states, rf, rp)
  s.events = events
  # Event order
  event_array = convert(Array{Float64, 2}, events)[:]
  event_order = sortperm(event_array)

  # Indicates that the first non -Inf event has been processed, and Δt can be calculated thereafter
  last_event_switch = false
  local last_event
  for i = 1:length(event_order)
    id, state_index = Tuple(CartesianIndices((events.individuals,
                                        length(_state_progressions[T][2:end])))[event_order[i]])
    new_state = _state_progressions[T][state_index+1]
    time = s.events[new_state][id]
    if time == -Inf
      @debug "Skipping event $i" Event{T}(time, id, new_state)
      continue
    elseif isnan(time)
      @debug "Loglikelihood calculation complete!"
      break
    end
    event = Event{T}(time, id, new_state)
    @debug "Calculating the loglikehood contribution of event $i" event
    # Get event rate totals
    rate_total = sum(s.event_rates)
    if rate_total < 0.0
      @error "Event rate total $i = $(round(rate_total, digits=3)), setting loglikelihood to -Inf"
      ll = -Inf
    elseif rate_total == 0.0
      @warn "Event rate total $i = $(round(rate_total, digits=3)), setting loglikelihood to -Inf"
      ll = -Inf
    else
      @debug "Event rate total $i = $(round(rate_total, digits=3))"
    end
    if loglikelihood_output
      if last_event_switch
        ΔT = event.time - last_event.time
        # Add event occurence contribution to loglikelihood
        ll += log(rate_total) - rate_total * ΔT
      end
      # Stop log likelihood calculation anytime the loglikelihood goes to -Inf
      if ll <= early_decision_value
        @debug "Loglikelihood calculation stopped early as loglikelihood <= decision value" loglikelihood = ll decision = early_decision_value
        if ll == -Inf
          @debug "Event $i resulted in a -Inf loglikelihood (transition of individual $id at t = $(round(time, digits=3)) to state $new_state)"
        end
        ll = -Inf
        break
      end
      # Get the individual rate associated with the event that occurred
      event_rate = s.event_rates[event.new_state][event.individual]
      # Add the specific event contribution to loglikelihood
      ll += log(event_rate/rate_total)
    end
    # Updates
    last_event = event
    if last_event_switch == false
      last_event_switch = true
    end
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
