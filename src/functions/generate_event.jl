"""
generate_event(rates::SEIR_Rates,
               time=0.::Float64)

A function to generate an event time and event type from a rate array
"""
function generate_event(rates::SEIR_Rates,
                        time=0.::Float64)
  totals = [sum(rates.exposure.external);
            sum(rates.exposure.internal);
            sum(rates.infection);
            sum(rates.removal)]
  total = sum(totals)
  if total == Inf
    time = time
    event_type = findfirst(totals .== Inf)
    event_index = findfirst(rates[event_type][:] .== Inf)
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


"""
generate_event(rates::SIR_Rates,
               time=0.::Float64)

A function to generate an event time and event type from a rate array
"""
function generate_event(rates::SIR_Rates,
                        time=0.::Float64)
  totals = [sum(rates.infection.external);
            sum(rates.infection.internal);
            sum(rates.removal)]
  total = sum(totals)
  if total == Inf
    time = time
    event_type = findfirst(totals .== Inf)
    event_index = findfirst(rates[event_type][:] .== Inf)
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
      event_index = findfirst(rand(Multinomial(1, rates.infection.external/totals[event_type])))
    elseif event_type == 2
      event_index = findfirst(rand(Multinomial(1, rates.infection.internal[:]/totals[event_type])))
    elseif event_type == 3
      event_index = findfirst(rand(Multinomial(1, rates.removal/totals[event_type])))
    end
  end
  return time, SIR_Event(event_type, event_index)
end


"""
generate_event(rates::SEI_Rates,
               time=0.::Float64)

A function to generate an event time and event type from a rate array
"""
function generate_event(rates::SEI_Rates,
                        time=0.::Float64)
  totals = [sum(rates.exposure.external);
            sum(rates.exposure.internal);
            sum(rates.infection)]
  total = sum(totals)
  if total == Inf
    time = time
    event_type = findfirst(totals .== Inf)
    event_index = findfirst(rates[event_type][:] .== Inf)
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
    end
  end
  return time, SEI_Event(event_type, event_index)
end


"""
generate_event(rates::SI_Rates,
               time=0.::Float64)

A function to generate an event time and event type from a rate array
"""
function generate_event(rates::SI_Rates,
                        time=0.::Float64)
  totals = [sum(rates.infection.external);
            sum(rates.infection.internal)]
  total = sum(totals)
  if total == Inf
    time = time
    event_type = findfirst(totals .== Inf)
    event_index = findfirst(rates[event_type][:] .== Inf)
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
      event_index = findfirst(rand(Multinomial(1, rates.infection.external/totals[event_type])))
    elseif event_type == 2
      event_index = findfirst(rand(Multinomial(1, rates.infection.internal[:]/totals[event_type])))
    end
  end
  return time, SI_Event(event_type, event_index)
end
