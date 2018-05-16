struct Simulation{T <: EpidemicModel}
  iterations::Int64
  population::DataFrame
  risk_functions::RiskFunctions{T}
  risk_parameters::RiskParameters{T}
  disease_states::Vector{DiseaseState}
  transmission_rates::TransmissionRates
  event_rates::EventRates{T}
  events::Events{T}
  transmission_network::TransmissionNetwork

  function Simulation{T}(pop::DataFrame,
                         rf::RiskFunctions{T},
                         rp::RiskParameters{T}) where T <: EpidemicModel

    n_ids = length(pop)

    # Initialize everything
    states = fill(State_S, n_ids)
    tr = initialize(TransmissionRates, states, pop, rf, rp)
    rates = initialize(EventRates{T}, tr, states, pop, rf, rp)
    events = Events{T}(n_ids)
    net = TransmissionNetwork(n_ids)

    # Generate initial event
    event = generate(Event{T}, rates)

    # Check to make sure initialization was valid
    event.time == Inf && error("Invalid simulation initialization")

    # Set time to 0.
    event.time = 0.

    # Update `Events` and `States`
    update!(events, event)
    update!(states, event)

    # If a transmission, generate the source before an update to `TransmissionRates` occurs
    if event.new_state in [State_E; State_I]
      xfer = generate(Transmission, tr)
      update!(net, xfer)
    end

    # Update `EventRates` before `TransmissionRates`
    update!(rates, tr, event, states, pop, rf, rp)
    update!(tr, event, states, pop, rf, rp)

    # Return initialized simulation objects
    return new{T}(1, pop, rf, rp, states, tr, rates, events, xfer)
  end
end
