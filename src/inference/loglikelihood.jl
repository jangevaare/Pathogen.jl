function loglikelihood(rp::RiskParameters{T},
                       rf::RiskFunctions{T},
                       events::Events{T},
                       pop::Population,
                       starting_states::DiseaseStates;
                       loglikelihood_output::Bool=true,
                       transmission_rates_output::Bool=true,
                       transmissions_output::Bool=true,
                       early_decision_value::Float64=-Inf,
                       transmission_rate_cache::Union{Nothing, TransmissionRateCache}=nothing,
                       workspace::Union{Nothing, LikelihoodWorkspace{T}}=nothing) where T <: DiseaseStateSequence
  n = individuals(pop)
  length(starting_states) == individuals(events) == n || throw(DimensionMismatch("Population, states and events must have equal lengths"))
  valid_states = convert(DiseaseStates, T)
  all(state -> state in valid_states, starting_states) || throw(ArgumentError("Invalid initial disease state"))
  work = workspace === nothing ? LikelihoodWorkspace(T, n) : workspace
  length(work.disease_states) == n || throw(DimensionMismatch("Likelihood workspace size differs from population"))
  copyto!(work.disease_states, starting_states)
  tr, rates = work.transmission_rates, work.event_rates
  _initialize!(tr, work.disease_states, pop, rf, rp)
  _initialize!(rates, tr, work.disease_states, pop, rf, rp)
  work.rate_total[] = sum(rates)
  empty!(work.transmissions)
  event_order = _event_order!(work, events)
  states = valid_states[2:end]
  ll = 0.0
  last_time = nothing
  for (i, index) in enumerate(event_order)
    id = mod1(index, n)
    new_state = states[fld(index-1, n)+1]
    time = events[new_state][id]
    time == -Inf && continue
    isnan(time) && break
    event = Event{T}(time, id, new_state)
    if _new_transmission(event) && transmissions_output
      push!(work.transmissions, id)
    end
    if loglikelihood_output
      # Refresh occasionally to bound drift from incremental total updates.
      if i % 256 == 0 || work.rate_total[] <= 0.0
        work.rate_total[] = sum(rates)
      end
      total = work.rate_total[]
      event_rate = rates[new_state][id]
      if !(total > 0.0 && event_rate > 0.0) || !isfinite(total) || !isfinite(event_rate)
        ll = -Inf
        break
      end
      if last_time === nothing
        ll += log(event_rate / total)
      else
        ll += log(event_rate) - total * (time - last_time)
      end
      # Partial continuous-time log densities are not acceptance bounds. The
      # threshold keyword remains accepted, but cannot safely prune this replay.
      if !isfinite(ll)
        ll = -Inf
        break
      end
    end
    last_time = time
    update!(work.disease_states, event)
    update!(tr, rates, event, work.disease_states, pop, rf, rp;
            transmission_rate_cache=transmission_rate_cache, rate_total=work.rate_total)
  end
  return loglikelihood_output ? ll : nothing,
         transmission_rates_output ? tr : nothing,
         transmissions_output ? work.transmissions : nothing
end
