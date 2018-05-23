mutable struct Simulation{T <: EpidemicModel}
  time::Float64
  iterations::Int64
  population::DataFrame
  risk_functions::RiskFunctions{T}
  risk_parameters::RiskParameters{T}
  disease_states::Vector{DiseaseState}
  transmission_rates::TransmissionRates
  event_rates::EventRates{T}
  events::Events{T}
  transmission_network::TransmissionNetwork

  function Simulation(pop::DataFrame,
                      rf::RiskFunctions{T},
                      rp::RiskParameters{T}) where T <: EpidemicModel

    n_ids = size(pop, 1)

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

    # Update `Events` and `States`
    update!(events, event)
    update!(states, event)

    # If a new transmission, generate the source before an update to `TransmissionRates` occurs
    if _new_transmission(event)
      update!(net, generate(Transmission,
                            tr,
                            event))
    end

    # Update `TransmissionRates` before `EventRates`
    update!(tr, event, states, pop, rf, rp)
    update!(rates, tr, event, states, pop, rf, rp)

    # Return initialized simulation objects
    return new{T}(0.0, 1, pop, rf, rp, states, tr, rates, events, net)
  end
end
