function initialize_simulation(pop::DataFrame,
                               rf::RiskFunctions{T},
                               rp::RiskParameters{T}) where T <: EpidemicModel

  n_ids = length(pop)

  # Initialize everything
  states = DiseaseStates{T}(n_ids)
  tr = initialize(::Type{TransmissionRates}, states, pop, rf, rp)
  rates = initialize(::Type{EventRates{T}}, tr, states, pop, rf, rp)
  events = Events{T}(n_ids)
  net = TransmissionNetwork(n_ids)

  # Generate initial event
  event = generate_event(rates)

  # Check to make sure initialization was valid
  event.time == Inf && error("Invalid simulation initialization")

  # Set time to 0.
  event.time = 0.

  # Update `Events` and `States`
  update!(events, event)
  update!(states, event)

  # If a transmission, generate the source before an update to `TransmissionRates` occurs
  if typeof(event) in [Exposure{SEIR}; Exposure{SEI}; Infection{SIR}; Infection{SI}]
    xfer = generate_transmission(tr)
    update!(net, xfer)
  end

  # Update `EventRates` before `TransmissionRates`
  update!(rates, tr, event, states, pop, rf, rp)
  update!(tr, event, states, pop, rf, rp)

  # Return initialized simulation objects
  return net, states, tr, rates, events
end
