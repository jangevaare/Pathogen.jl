function generate(::Type{Event{T}},
                  rates::EventRates{T},
                  time=0.::Float64) where T <: EpidemicModel

  if T == SEIR
    totals = [sum(rates.exposure); sum(rates.infection); sum(rates.removal)]
  elseif T == SEI
    totals = [sum(rates.exposure); sum(rates.infection); 0.]
  elseif T == SIR
    totals = [0.; sum(rates.infection); sum(rates.removal)]
  elseif T == SI
    totals = [0.; sum(rates.infection); 0.]
  end
  cumtotals = cumsum(totals)
  if cumtotals[end] == Inf
    time = time
    present = [State_I; State_E; State_R][findfirst(cumtotals .== Inf)]
    id = findfirst(rates[present] .== Inf)
  elseif total == 0.
    time = Inf
    event_type = 0
    event_index = 0
  else
    # Generate event time
    time += rand(Exponential(1/total))
    # Generate event type and index
    event_type = findfirst(rand(Multinomial(1, totals/total)))
    if event_type == 1
      event_index = findfirst(rand(Multinomial(1, rates.exposure.external/totals[event_type])))
    elseif event_type == 2
      event_index = findfirst(rand(Multinomial(1, rates.exposure.internal[:]/totals[event_type])))
    elseif event_type == 3
      event_index = findfirst(rand(Multinomial(1, rates.infection/totals[event_type])))
    elseif event_type == 4
      event_index = findfirst(rand(Multinomial(1, rates.removal/totals[event_type])))
    end
  end
  return time, SEIR_Event(event_type, event_index)
end
