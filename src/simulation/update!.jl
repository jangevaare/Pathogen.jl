function update!(sim::Simulation{T, M},
                 event::Event{T},
                 transmission::Transmission) where {T <: DiseaseStateSequence, M <: ILM}
  if _time(event) < sim.time
    error("Time of a new event must be >= that of the previous event")
  elseif _time(event) < Inf
    update!(sim.events, event)
    update!(sim.transmission_network, transmission)
    update!(sim.disease_states, event)
    update!(sim.transmission_rates,
            sim.event_rates,
            event,
            sim.disease_states,
            sim.population,
            sim.risk_functions,
            sim.risk_parameters)
    # Count the iteration as having occurred
    sim.iterations += 1
  end
  # Update simulation time
  sim.time = _time(event)
  return sim
end