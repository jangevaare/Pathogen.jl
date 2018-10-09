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
  # To cause event ordering to ignore event times prior to starting_time (including -Inf events that led to starting states) set to NaN
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
    rate_total = sum(sum(s.event_rates[new_state]) for state in _state_progressions[T][2:end])
    if rate_total <= 0.0
      @error "Event rate total $i = $(round(rate_total, digits=3)) (it should be > 0.0)"
      ll = -Inf
    else
      @debug "Event rate total $i = $(round(rate_total, digits=3))"
    end
    if loglikelihood_output
      if last_event_switch
        ΔT = event.time - last_event.time
        # Add event occurence contribution to loglikelihood
        ll += log(rate_total) - rate_total * ΔT
      else
        last_event_switch = true
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
                       loglikelihood_output::Bool=true,
                       transmission_network_output::Bool=true,
                       early_decision_value::Float64=-Inf) where T <: EpidemicModel
  return loglikelihood(rp, rf, events, pop,
                       fill(State_S, pop.individuals),
                       loglikelihood_output = loglikelihood_output,
                       transmission_network_output = transmission_network_output,
                       early_decision_value = early_decision_value)
end


function loglikelihood2(rp::RiskParameters{T},
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
  # To cause event ordering to ignore event times prior to starting_time (including -Inf events that led to starting states) set to NaN
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
      break
    end
    event = Event{T}(time, id, new_state)
    @debug "Calculating the loglikehood contribution of event $i" event
    # Get exposure or infection event rate totals
    rate_total = sum(s.transmission_rates)
    if rate_total < 0.0
      @error "Rate total must be positive" rate_total
    elseif rate_total == 0.0
      @debug "Rate total is 0.0, skip ahead"
    else
      @debug "Event rate total $i = $(round(rate_total, digits=3))"
      if loglikelihood_output
        if last_event_switch
          ΔT = event.time - last_event.time
            # Add event occurence contribution to loglikelihood
          if _new_transmission(event)
            # log of exponential distribution's PDF
            ll += log(rate_total) - rate_total * ΔT
          else
            # log of exponential distribution's CCDF
            ll += -rate_total * ΔT
          end
        else
          last_event_switch = true
        end
        # Stop loglikelihood calculation anytime the loglikelihood goes to -Inf
        if ll <= early_decision_value
          @debug "Loglikelihood calculation stopped early as loglikelihood <= decision value" loglikelihood = ll decision = early_decision_value
          if ll == -Inf
            @debug "Event $i resulted in a -Inf loglikelihood (transition of individual $id at t = $(round(time, digits=3)) to state $new_state)"
          end
          ll = -Inf
          break
        end
        # Get the individual rate associated with the event that occurred
        event_rate = sum(s.transmission_rates, event.individual)
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
            event,
            s.disease_states,
            s.population,
            s.risk_functions,
            s.risk_parameters)
  end
  @debug "Loglikelihood calculation complete!"
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
