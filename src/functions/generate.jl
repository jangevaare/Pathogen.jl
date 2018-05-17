function generate(::Type{Event{T}},
                  rates::EventRates{T},
                  time=0.::Float64) where T <: EpidemicModel
  totals = [sum(rates[state]) for state in _state_progressions[T][2:end]]
  ctotals = cumsum(totals)
  if ctotals[end] == Inf
    time = time
    new_state = _state_progressions[T][findfirst(ctotals .== Inf) + 1]
    id = findfirst(rates[new_state] .== Inf)
  elseif ctotals[end] == 0.
    time = Inf
  else
    # Generate event time
    time += rand(Exponential(1/ctotals[end]))
    # Generate new state
    new_state_index = findfirst(rand(Multinomial(1, totals ./ ctotals[end])))
    new_state = _state_progressions[T][new_state_index+1]
    # Generate event indvidual
    id = findfirst(rand(Multinomial(1, rates[new_state] ./ totals[new_state_index])))
    return Event{T}(time, id, new_state)
  end
end

function generate(::Type{Transmission}, tr::TransmissionRates, id::Int64)
  external_or_internal = [tr.external[id]; sum(tr.internal[:,id])]
  if rand() <= external_or_internal[1]/sum(external_or_internal)
    return ExogenousTransmission(id)
  else
    source = findfirst(rand(Multinomial(1, tr.internal[:,id] ./ external_or_internal[2])))
    return EndogenousTransmission(id, source)
  end
end
