mutable struct Simulation{S <: DiseaseStateSequence, M <: ILM}
  time::Float64
  iterations::Int64
  population::Population
  risk_functions::RiskFunctions{S}
  risk_parameters::RiskParameters{S}
  disease_states::DiseaseStates
  transmission_rates::TransmissionRates
  event_rates::EventRates{S}
  events::Events{S}
  transmission_network::TransmissionNetwork

  function Simulation{M}(pop::Population,
                         rf::RiskFunctions{S},
                         rp::RiskParameters{S}) where {
                         S <: DiseaseStateSequence, 
                         M <: ILM}
    states = fill(State_S, pop.individuals)
    tr = initialize(TransmissionRates, states, pop, rf, rp)
    rates = initialize(EventRates, tr, states, pop, rf, rp)
    events = Events{S}(pop.individuals)
    net = TransmissionNetwork(pop.individuals)
    return new{S, M}(0.0, 0, pop, rf, rp, states, tr, rates, events, net)
  end


  function Simulation{M}(pop::Population,
                        states::DiseaseStates,
                        time::Float64,
                        rf::RiskFunctions{S},
                        rp::RiskParameters{S}) where {
                        S <: DiseaseStateSequence, 
                        M <: ILM}
    @debug "Initializing $T Simulation with the following starting states:" states
    if length(states) != pop.individuals
      @error "Length of initial disease state vector must match number of individuals"
    elseif !all(in.(states, Ref(convert(DiseaseStates, S))))
      @error "All states in initial disease state vector must be valid within specified epidemic model"
    end
    tr = initialize(TransmissionRates, states, pop, rf, rp)
    rates = initialize(EventRates, tr, states, pop, rf, rp)
    events = Events{S}(states)
    net = TransmissionNetwork(pop.individuals)
    return new{S, M}(time, 0, pop, rf, rp, copy(states), tr, rates, events, net)
  end

  function Simulation{M}(pop::Population,
                         states::DiseaseStates,
                         rf::RiskFunctions{S},
                         rp::RiskParameters{S}) where {
                         S <: DiseaseStateSequence, 
                         M <: ILM}
    @debug "Initializing $T Simulation with the following starting states:" states
    if length(states) != pop.individuals
      @error "Length of initial disease state vector must match number of individuals"
    elseif !all(in.(states, Ref(convert(DiseaseStates, S))))
      @error "All states in initial disease state vector must be valid within specified epidemic model"
    end
    tr = initialize(TransmissionRates, states, pop, rf, rp)
    rates = initialize(EventRates, tr, states, pop, rf, rp)
    events = Events{S}(states)
    net = TransmissionNetwork(pop.individuals)
    return new{S, M}(0.0, 0, pop, rf, rp, copy(states), tr, rates, events, net)
  end
end

function Base.show(io::IO, 
                   x::Simulation{S, M}) where {
                   S <: DiseaseStateSequence, 
                   M <: ILM}
  y = "$S $M epidemic simulation @ time = $(round(x.time, digits=2))\n"

  
  for s in convert(DiseaseStates, S)
    y *= "\n$s = $(sum(x.disease_states .== Ref(s)))"
  end
  return print(io, y)
end