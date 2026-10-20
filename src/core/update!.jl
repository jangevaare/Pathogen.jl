function update!(events::Events{T},
                 event::Event{T}) where T <: DiseaseStateSequence
  events[event.new_state][event.individual] = event.time
  return events
end

function update!(states::DiseaseStates,
                 event::Event{T}) where T <: DiseaseStateSequence
  states[event.individual] = event.new_state
  return states
end

function update!(net::TransmissionNetwork,
                 xfer::ExogenousTransmission)
  net.external[xfer.individual] = true
  return net
end

function update!(net::TransmissionNetwork,
                 xfer::EndogenousTransmission)
  net.internal[xfer.source, xfer.individual] = true
  return net
end

function update!(net::TransmissionNetwork,
                 xfer::NoTransmission)
  return net
end

function update!(tr::TransmissionRates,
                 event::Event{T}, states::DiseaseStates, pop::Population,
                 rf::RiskFunctions{T}, rp::RiskParameters{T};
                 xzero::Bool=false,
                 transmission_rate_cache::Union{Nothing, TransmissionRateCache}=nothing) where T <: DiseaseStateSequence
  id = event.individual
  if _new_transmission(event) && xzero
    tr.external[id] = 0.0
    tr.internal[:, id] .= 0.0
  end
  if event.new_state == State_I
    transmissibility = nothing
    for i in eachindex(states)
      states[i] == State_S || continue
      cached = transmission_rate_cache === nothing ? nothing : transmission_rate_cache.internal[id, i]
      if cached === nothing
        if transmissibility === nothing
          transmissibility = rf.transmissibility(rp.transmissibility, pop, id)
        end
        tr.internal[id, i] = rf.susceptibility(rp.susceptibility, pop, i) *
                             rf.infectivity(rp.infectivity, pop, i, id) * transmissibility
        if transmission_rate_cache !== nothing
          transmission_rate_cache.internal[id, i] = tr.internal[id, i]
        end
      else
        tr.internal[id, i] = cached
      end
    end
  elseif event.new_state == State_R
    for i in eachindex(states)
      if states[i] == State_S
        tr.internal[id, i] = 0.0
      end
    end
  end
  return tr
end

# Return the change in the global hazard while updating individual hazards.
# Standalone calls retain the historical convention that tr is updated first.
function _update_event_rates!(rates::EventRates{T}, tr::TransmissionRates,
                              event::Event{T}, states::DiseaseStates,
                              pop::Population, rf::RiskFunctions{T},
                              rp::RiskParameters{T}) where T <: DiseaseStateSequence
  id = event.individual
  delta = 0.0
  if event.new_state == State_E
    delta -= rates.exposure[id]
    rates.exposure[id] = 0.0
    rates.infection[id] = rf.latency(rp.latency, pop, id)
    delta += rates.infection[id]
  elseif event.new_state == State_I
    delta -= rates.infection[id]
    rates.infection[id] = 0.0
    if State_R ∈ T
      rates.removal[id] = rf.removal(rp.removal, pop, id)
      delta += rates.removal[id]
    end
    target = State_E ∈ T ? rates.exposure : rates.infection
    for i in eachindex(states)
      if states[i] == State_S
        target[i] += tr.internal[id, i]
        delta += tr.internal[id, i]
      end
    end
  elseif event.new_state == State_R
    delta -= rates.removal[id]
    rates.removal[id] = 0.0
    target = State_E ∈ T ? rates.exposure : rates.infection
    for i in eachindex(states)
      if states[i] == State_S
        updated = sum(tr, i)
        delta += updated - target[i]
        target[i] = updated
      end
    end
  end
  return delta
end

function update!(rates::EventRates{T}, tr::TransmissionRates,
                 event::Event{T}, states::DiseaseStates, pop::Population,
                 rf::RiskFunctions{T}, rp::RiskParameters{T}) where T <: DiseaseStateSequence
  _update_event_rates!(rates, tr, event, states, pop, rf, rp)
  return rates
end

function update!(tr::TransmissionRates, rates::EventRates{T},
                 event::Event{T}, states::DiseaseStates, pop::Population,
                 rf::RiskFunctions{T}, rp::RiskParameters{T};
                 xzero::Bool=false,
                 transmission_rate_cache::Union{Nothing, TransmissionRateCache}=nothing,
                 rate_total::Union{Nothing, Base.RefValue{Float64}}=nothing) where T <: DiseaseStateSequence
  if event.new_state == State_R
    # Subtract before clearing; preserve historical columns for infected people.
    target = State_E ∈ T ? rates.exposure : rates.infection
    delta = -rates.removal[event.individual]
    no_sources = !any(==(State_I), states)
    for i in eachindex(states)
      if states[i] == State_S
        previous = target[i]
        departing = tr.internal[event.individual, i]
        tr.internal[event.individual, i] = 0.0
        # Reset exactly when the final source leaves, avoiding residual hazard.
        updated = no_sources ? tr.external[i] : max(tr.external[i], previous - departing)
        if !no_sources && departing > 0.0 && updated <= sqrt(eps(Float64)) * previous
          # Subtracting a dominant source can erase a much smaller surviving
          # source. Rebuild only the affected column when cancellation is severe.
          updated = sum(tr, i)
        end
        delta += updated - previous
        target[i] = updated
      end
    end
    rates.removal[event.individual] = 0.0
  else
    update!(tr, event, states, pop, rf, rp; xzero=xzero, transmission_rate_cache=transmission_rate_cache)
    delta = _update_event_rates!(rates, tr, event, states, pop, rf, rp)
  end
  if rate_total !== nothing
    previous = rate_total[]
    updated = previous + delta
    rate_total[] = delta < 0.0 && updated <= sqrt(eps(Float64)) * previous ? sum(rates) : updated
  end
  return tr, rates
end
