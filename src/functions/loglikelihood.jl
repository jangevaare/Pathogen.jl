function loglikelihood(rp::RiskParameters{T},
                       rf::RiskFunctions{T},
                       events::Events{T},
                       pop::Population,
                       starting_states::Vector{DiseaseState};
                       starting_time::Float64=0.0,
                       loglikelihood_output::Bool=true,
                       transmission_network_output::Bool=true,
                       early_decision_value::Float64=-Inf) where T <: EpidemicModel
  # Initialize
  ll = 0.
  # We can use the Simulation struct to recreate epidemic
  s = Simulation(pop, starting_states, starting_time, rf, rp)
  s.events = events
  # Event order
  event_array = convert(Array{Float64, 2}, events)[:]
  # To cause event ordering to ignore -Inf event times (i.e. events that led to starting states) set to NaN
  event_array[event_array .== -Inf] .= NaN
  event_order = sortperm(event_array)
  local last_event
  local ΔT
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
      if i == 1
        ΔT = event.time - starting_time
      else
        ΔT = event.time - last_event.time
      end
      # Add event occurence contribution to loglikelihood
      ll += log(rate_total) - rate_total * ΔT
      # Stop log likelihood calculation anytime the loglikelihood goes to -Inf
      if ll <= early_decision_value
        if ll == -Inf
          @debug "Event $i resulted in a -Inf loglikelihood (transition of individual $id at t = $(round(time, 3)) to state $new_state)"
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

function loglikelihood(rp::RiskParameters{T},
                       rf::RiskFunctions{T},
                       events::Events{T},
                       pop::Population;
                       starting_time::Float64=0.0,
                       loglikelihood_output::Bool=true,
                       transmission_network_output::Bool=true,
                       early_decision_value::Float64=-Inf) where T <: EpidemicModel
  return loglikelihood(rp, rf, events, pop,
                       fill(State_S, pop.individuals),
                       starting_time = starting_time,
                       loglikelihood_output = loglikelihood_output,
                       transmission_network_output = transmission_network_output,
                       early_decision_value = early_decision_value)
end
