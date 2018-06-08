function next!(s::Simulation{T}) where T <: EpidemicModel
  event = generate(Event{T}, s.event_rates, s.time)
  if event.time < Inf
    push!(s.events, event)
    update!(s.disease_states, event)
    # Checks if a _new_transmission within function
    # Could make this into a single function for conveinence (also appears in Simulation.jl)
    update!(s.transmission_network,
            generate(Transmission,
                     s.transmission_rates,
                     event))
    # Update `EventRates` before `TransmissionRates`
    update!(s.event_rates,
            s.transmission_rates,
            event,
            s.disease_states,
            s.population,
            s.risk_functions,
            s.risk_parameters)
    update!(s.transmission_rates,
            event,
            s.disease_states,
            s.population,
            s.risk_functions,
            s.risk_parameters)
    # Count the iteration as having occurred.
    s.iterations += 1
  end
  # Update simulation time
  s.time = event.time
  # Return updated Simulation object
  return s
end
