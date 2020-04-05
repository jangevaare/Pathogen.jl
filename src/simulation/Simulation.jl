mutable struct Simulation{S <: DiseaseStateSequence, M <: ILM}
  start_time::Float64
  simulation_time::Float64
  iterations::Int64
  population::Population
  risk_functions::RiskFunctions{S}
  risk_parameters::RiskParameters{S}
  substitution_model::Union{Nothing, NucleicAcidSubstitutionModel}
  disease_states::DiseaseStates
  transmission_rates::TransmissionRates
  event_rates::EventRates{S}
  events::Events{S}
  transmission_network::TransmissionNetwork

  function Simulation{S, M}(pop::Population,
                            rf::RiskFunctions{S},
                            rp::RiskParameters{S}) where {
                            S <: DiseaseStateSequence, 
                            M <: TNILM}
    states = fill(State_S, pop.individuals)
    tr = initialize(TransmissionRates, states, pop, rf, rp)
    rates = initialize(EventRates, tr, states, pop, rf, rp)
    events = Events{S}(pop.individuals)
    net = TransmissionNetwork(states)
    return new{S, M}(0.0, 0.0, 0, pop, rf, rp, nothing, states, tr, rates, events, net)
  end


  function Simulation{S, M}(pop::Population,
                            states::DiseaseStates,
                            time::Float64,
                            rf::RiskFunctions{S},
                            rp::RiskParameters{S}) where {
                            S <: DiseaseStateSequence, 
                            M <: TNILM}
    @debug "Initializing $S $M Simulation with the following starting states:" states
    if length(states) != pop.individuals
      throw(ErrorException("Length of initial disease state vector must match number of individuals"))
    elseif !all(in.(states, Ref(convert(DiseaseStates, S))))
      throw(ErrorException("All states in initial disease state vector must be valid within specified epidemic model"))
    end
    tr = initialize(TransmissionRates, states, pop, rf, rp)
    rates = initialize(EventRates, tr, states, pop, rf, rp)
    events = Events{S}(states)
    net = TransmissionNetwork(states)
    return new{S, M}(time, time, 0, pop, rf, rp, nothing, copy(states), tr, rates, events, net)
  end

  function Simulation{S, M}(pop::Population,
                            states::DiseaseStates,
                            rf::RiskFunctions{S},
                            rp::RiskParameters{S}) where {
                            S <: DiseaseStateSequence, 
                            M <: TNILM}
    @debug "Initializing $S $M Simulation with the following starting states:" states
    if length(states) != pop.individuals
      throw(ErrorException("Length of initial disease state vector must match number of individuals"))
    elseif !all(in.(states, Ref(convert(DiseaseStates, S))))
      throw(ErrorException("All states in initial disease state vector must be valid within specified epidemic model"))
    end
    tr = initialize(TransmissionRates, states, pop, rf, rp)
    rates = initialize(EventRates, tr, states, pop, rf, rp)
    events = Events{S}(states)
    net = TransmissionNetwork(states)
    return new{S, M}(0.0, 0.0, 0, pop, rf, rp, nothing, copy(states), tr, rates, events, net)
  end

function Simulation{S, M}(pop::Population,
                          rf::RiskFunctions{S},
                          rp::RiskParameters{S},
                          sm::NASM) where {
                          S <: DiseaseStateSequence, 
                          M <: PhyloILM}
    states = fill(State_S, pop.individuals)
    tr = initialize(TransmissionRates, states, pop, rf, rp)
    rates = initialize(EventRates, tr, states, pop, rf, rp)
    events = Events{S}(pop.individuals)
    net = TransmissionNetwork(states)
    return new{S, M}(0.0, 0.0, 0, pop, rf, rp, sm, states, tr, rates, events, net)
  end


  function Simulation{S, M}(pop::Population,
                            states::DiseaseStates,
                            time::Float64,
                            rf::RiskFunctions{S},
                            rp::RiskParameters{S},
                            sm::NASM) where {
                            S <: DiseaseStateSequence, 
                            M <: PhyloILM}
    @debug "Initializing $S $M Simulation with the following starting states:" states
    if length(states) != pop.individuals
      throw(ErrorException("Length of initial disease state vector must match number of individuals"))
    elseif !all(in.(states, Ref(convert(DiseaseStates, S))))
      throw(ErrorException("All states in initial disease state vector must be valid within specified epidemic model"))
    end
    tr = initialize(TransmissionRates, states, pop, rf, rp)
    rates = initialize(EventRates, tr, states, pop, rf, rp)
    events = Events{S}(states)
    net = TransmissionNetwork(states)
    return new{S, M}(time, time, 0, pop, rf, rp, sm, copy(states), tr, rates, events, net)
  end

  function Simulation{S, M}(pop::Population,
                            states::DiseaseStates,
                            rf::RiskFunctions{S},
                            rp::RiskParameters{S},
                            sm::NASM) where {
                            S <: DiseaseStateSequence, 
                            M <: PhyloILM}
    @debug "Initializing $S $M Simulation with the following starting states:" states
    if length(states) != pop.individuals
      throw(ErrorException("Length of initial disease state vector must match number of individuals"))
    elseif !all(in.(states, Ref(convert(DiseaseStates, S))))
      throw(ErrorException("All states in initial disease state vector must be valid within specified epidemic model"))
    end
    tr = initialize(TransmissionRates, states, pop, rf, rp)
    rates = initialize(EventRates, tr, states, pop, rf, rp)
    events = Events{S}(states)
    net = TransmissionNetwork(states)
    return new{S, M}(0.0, 0.0, 0, pop, rf, rp, sm, copy(states), tr, rates, events, net)
  end
end

function Base.show(io::IO, 
                   x::Simulation{S, M}) where {
                   S <: DiseaseStateSequence, 
                   M <: ILM}
  y = "$S $M epidemic simulation @ time = $(round(x.simulation_time, digits=2))\n"

  for s in convert(DiseaseStates, S)
    y *= "\n$s = $(sum(x.disease_states .== Ref(s)))"
  end
  return print(io, y)
end