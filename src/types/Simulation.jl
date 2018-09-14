mutable struct Simulation{T <: EpidemicModel}
  time::Float64
  iterations::Int64
  population::Population
  risk_functions::RiskFunctions{T}
  risk_parameters::RiskParameters{T}
  disease_states::Vector{DiseaseState}
  transmission_rates::TransmissionRates
  event_rates::EventRates{T}
  events::Events{T}
  transmission_network::TransmissionNetwork

  function Simulation(pop::Population,
                      rf::RiskFunctions{T},
                      rp::RiskParameters{T}) where T <: EpidemicModel
    states = fill(State_S, pop.individuals)
    tr = initialize(TransmissionRates, states, pop, rf, rp)
    rates = initialize(EventRates{T}, tr, states, pop, rf, rp)
    events = Events{T}(pop.individuals)
    net = TransmissionNetwork(pop.individuals)
    return new{T}(0.0, 0, pop, rf, rp, states, tr, rates, events, net)
  end

  function Simulation(pop::Population,
                      states::Vector{DiseaseState},
                      time::Float64,
                      rf::RiskFunctions{T},
                      rp::RiskParameters{T}) where T <: EpidemicModel
    if length(states) != pop.individuals
      @error "Length of initial disease state vector must match number of individuals"
    elseif !all(in.(states, Ref(_state_progressions[T])))
      @error "All states in initial disease state vector must be valid within specified epidemic model"
    end
    tr = initialize(TransmissionRates, states, pop, rf, rp)
    rates = initialize(EventRates{T}, tr, states, pop, rf, rp)
    events = Events{T}(states)
    net = TransmissionNetwork(pop.individuals)
    return new{T}(time, 0, pop, rf, rp, states, tr, rates, events, net)
  end
end

function Base.show(io::IO, x::Simulation{T}) where T <: EpidemicModel
  return print(io, "$T epidemic simulation")
end
