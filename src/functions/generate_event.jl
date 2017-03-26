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
    event_index = findfirst(rand(Multinomial(1, rates[event_type][:]/totals[event_type])))
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
    event_index = findfirst(rand(Multinomial(1, rates[event_type][:]/totals[event_type])))
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
    event_index = findfirst(rand(Multinomial(1, rates[event_type][:]/totals[event_type])))
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
    event_index = findfirst(rand(Multinomial(1, rates[event_type][:]/totals[event_type])))
  end
  return time, SI_Event(event_type, event_index)
end
