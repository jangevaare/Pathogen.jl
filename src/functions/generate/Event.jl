function generate(::Type{Event},
                  rates::EventRates{T},
                  time::Float64) where T <: EpidemicModel
  totals = Weights([sum(rates[state]) for state in _state_progressions[T][2:end]])
  if sum(totals) == Inf
    new_state = sample(_state_progressions[T][2:end], totals)
    id = sample(1:rates.individuals, Weights(rates[new_state]))
    return Event{T}(time, id, new_state)
  elseif sum(totals) == 0.0
    return Event{T}(Inf)
  else
    # Generate new state
    new_state = sample(_state_progressions[T][2:end], totals)
    # Generate event individual
    id = sample(1:rates.individuals, Weights(rates[new_state]))
    return Event{T}(time + rand(Exponential(1.0/sum(totals))), id, new_state)
  end
end

function generate(::Type{Event},
                  last_event::Event{T},
                  σ::Float64,
                  obs::EventObservations{T},
                  extents::EventExtents{T},
                  events::Events{T}) where T <: EpidemicModel
  # Calculate the event bounds without any transmission network restrictions
  lb, ub = _bounds(last_event, extents, obs, events)
  time = rand(TruncatedNormal(last_event.time,
                              σ,
                              lb,
                              ub))
  return Event(time, last_event)
end

function generate(::Type{Event},
                  last_event::Event{T},
                  σ::Float64,
                  obs::EventObservations{T},
                  extents::EventExtents{T},
                  events::Events{T},
                  network::TransmissionNetwork) where T <: EpidemicModel
  # Calculate the event bounds including transmission network restrictions
  lb, ub = _bounds(last_event, extents, obs, events, network)
  time = rand(TruncatedNormal(last_event.time,
                              σ,
                              lb,
                              ub))
  return Event(time, last_event)
end

function generate(::Type{Event},
                  last_event::Event{T},
                  rp::RiskParameters{T},
                  rf::RiskFunctions{T},
                  pop::Population,
                  obs::EventObservations{T},
                  extents::EventExtents{T},
                  events::Events{T}) where T <: EpidemicModel
  # Gibbs-type new `Event` generation
  # Calculate the event bounds without any transmission network restrictions
  lb, ub = _bounds(last_event, extents, obs, events)
  if (last_event.new_state == State_I) && (T in [SEIR; SEI])
    λ = rf.latency(rp.latency, pop, last_event.individual)
    time = rand(Truncated(Exponential(1.0/λ), lb, ub))
  elseif last_event.new_state == State_R
    λ = rf.removal(rp.removal, pop, last_event.individual)
    time = rand(Truncated(Exponential(1.0/λ), lb, ub))
  else
    @error "Gibbs event time generation is only supported for latent and infectious periods"
  end
  return Event(time, last_event)
end

function generate(::Type{Event},
                  last_event::Event{T},
                  rp::RiskParameters{T},
                  rf::RiskFunctions{T},
                  pop::Population,
                  obs::EventObservations{T},
                  extents::EventExtents{T},
                  network::TransmissionNetwork) where T <: EpidemicModel
  # Gibbs-type new `Event` generation
  # Calculate the event bounds including transmission network restrictions
  lb, ub = _bounds(last_event, extents, obs, events, network)
  if (last_event.new_state == State_I) && (T in [SEIR; SEI])
    λ = rf.latency(rp.latency, pop, last_event.individual)
    time = rand(Truncated(Exponential(1.0/λ), lb, ub))
  elseif last_event.new_state == State_R
    λ = rf.removal(rp.removal, pop, last_event.individual)
    time = rand(Truncated(Exponential(1.0/λ), lb, ub))
  else
    @error "Gibbs event time generation is only supported for latent and infectious periods"
  end
  return Event(time, last_event)
end
