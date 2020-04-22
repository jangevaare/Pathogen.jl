function generate(::Type{Event},
                  rates::EventRates{T},
                  time::Float64) where T <: EpidemicModel
  totals = Weights([sum(rates[state]) for state in _state_progressions[T][2:end]])
  if sum(totals) == Inf
    new_state = sample(_state_progressions[T][2:end], totals)
    id = sample(1:rates.individuals, Weights(rates[new_state]))
    return Event{T}(time, id, new_state)
  elseif sum(totals) == 0.0
    return NoEvent()
  else
    # Generate new state
    new_state = sample(_state_progressions[T][2:end], totals)
    # Generate event individual
    id = sample(1:rates.individuals, Weights(rates[new_state]))
    return Event{T}(time + rand(Exponential(1.0 / sum(totals))), id, new_state)
  end
end

function generate(::Type{Transmission},
                  tr::TransmissionRates,
                  event::Event{T}) where T <: EpidemicModel
  id = event.individual
  if _new_transmission(event)
    external_or_internal = Weights([tr.external[id]; sum(tr.internal[:,id])])
    if sum(external_or_internal) == 0.0
      @error "All transmission rates = 0.0, No transmission can be generated"
      return NoTransmission()
    elseif sample([true; false], external_or_internal)
      @debug "Exogenous tranmission generated"
      return ExogenousTransmission(id)
    else
      source = sample(1:individuals(tr), Weights(tr.internal[:, id]))
      @debug "Endogenous transmission generated (source id = $source)"
      return EndogenousTransmission(id, source)
    end
  else
    @debug "No transmission generated"
    return NoTransmission()
  end
end