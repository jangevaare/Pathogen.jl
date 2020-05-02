function loglikelihood(rp::RiskParameters{S},
                       rf::RiskFunctions{S},
                       events::Events{S},
                       pop::Population,
                       starting_states::DiseaseStates;
                       loglikelihood_output::Bool=true,
                       transmission_rates_output::Bool=true,
                       transmissions_output::Bool=true,
                       early_decision_value::Float64=-Inf) where {
                       S <: DiseaseStateSequence}
  # Initialize
  ll = 0.0
  transmissions = Int64[]
  # We can use the Simulation struct to recreate epidemic
  s = Simulation{S, TNILM}(pop, starting_states, rf, rp)
  s.events = events
  # Event order
  event_array = convert(Array{Float64, 2}, events)[:]
  event_order = sortperm(event_array)

  # Indicates that the first non -Inf event has been processed, and Δt can be calculated thereafter
  last_event_switch = false
  local last_event
  for i = 1:length(event_order)
    id, state_index = Tuple(CartesianIndices((individuals(events),
                                        length(convert(DiseaseStates, S)[2:end])))[event_order[i]])
    new_state = convert(DiseaseStates, S)[state_index+1]
    time = s.events[new_state][id]
    if time == -Inf
      @debug "Skipping initial event $i" Event{S}(time, id, new_state)
      continue
    elseif isnan(time)
      @debug "Loglikelihood calculation complete!" loglikelihood = ll
      break
    end
    event = Event{S}(time, id, new_state)
    if _new_transmission(event) & transmissions_output
      push!(transmissions, id)
    end
    @debug "Calculating the loglikehood contribution of event $i" event
    # Get event rate totals
    rate_total = sum(s.event_rates)
    if rate_total < 0.0
      @error "Event rate total $i = $(round(rate_total, digits=3)), setting loglikelihood to -Inf"
      ll = -Inf
    elseif rate_total == 0.0
      @debug "Event rate total $i = $(round(rate_total, digits=3)), setting loglikelihood to -Inf"
      ll = -Inf
    else
      @debug "Event rate total $i = $(round(rate_total, digits=3))"
    end
    if loglikelihood_output
      if last_event_switch
        ΔT = _time(event) - _time(last_event)
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
      ll += log(event_rate / rate_total)
    end
    # Updates
    last_event = event
    if last_event_switch == false
      last_event_switch = true
    end
    update!(s.disease_states, event)
    update!(s.transmission_rates,
            s.event_rates,
            event,
            s.disease_states,
            s.population,
            s.risk_functions,
            s.risk_parameters)
  end
  return loglikelihood_output ? ll : nothing,
         transmission_rates_output ? s.transmission_rates : nothing,
         transmissions_output ? transmissions : nothing
end
