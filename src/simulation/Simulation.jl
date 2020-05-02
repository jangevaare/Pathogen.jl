mutable struct Simulation{S <: DiseaseStateSequence, M <: ILM}
  time::Float64
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
                            rp::RiskParameters{S};
                            time::Float64 = 0.0) where {
                            S <: DiseaseStateSequence,
                            M <: TNILM}
    states = fill(State_S, individuals(pop))
    tr = initialize(TransmissionRates, states, pop, rf, rp)
    rates = initialize(EventRates, tr, states, pop, rf, rp)
    events = Events{S}(individuals(pop))
    net = TransmissionNetwork(states)
    return new{S, M}(time, 0, pop, rf, rp, nothing, states, tr, rates, events, net)
  end

  function Simulation{S, M}(pop::Population,
                            states::DiseaseStates,
                            rf::RiskFunctions{S},
                            rp::RiskParameters{S};
                            time::Float64 = 0.0,
                            skip_checks::Bool=false) where {
                            S <: DiseaseStateSequence,
                            M <: TNILM}
    @debug "Initializing $S $M Simulation with the following starting states:" states = convert(DiseaseStates, S) counts = [sum(states .== s) for s in convert(DiseaseStates, S)]
    if !skip_checks
      if length(states) != individuals(pop)
        @error "Length of initial disease state vector must match number of individuals"
      elseif !all(in.(states, Ref(convert(DiseaseStates, S))))
        @error "All states in initial disease state vector must be valid within specified epidemic model"
      end
    end
    tr = initialize(TransmissionRates, states, pop, rf, rp)
    rates = initialize(EventRates, tr, states, pop, rf, rp)
    events = Events{S}(states)
    net = TransmissionNetwork(states)
    return new{S, M}(time, 0, pop, rf, rp, nothing, copy(states), tr, rates, events, net)
  end

  function Simulation{S, M}(pop::Population,
                            rf::RiskFunctions{S},
                            rp::RiskParameters{S},
                            sm::NASM;
                            time::Float64 = 0.0) where {
                            S <: DiseaseStateSequence,
                            M <: PhyloILM}
    states = fill(State_S, individuals(pop))
    tr = initialize(TransmissionRates, states, pop, rf, rp)
    rates = initialize(EventRates, tr, states, pop, rf, rp)
    events = Events{S}(individuals(pop))
    net = TransmissionNetwork(states)
    return new{S, M}(time, 0, pop, rf, rp, sm, states, tr, rates, events, net)
  end

  function Simulation{S, M}(pop::Population,
                            states::DiseaseStates,
                            rf::RiskFunctions{S},
                            rp::RiskParameters{S},
                            sm::NASM;
                            time::Float64 = 0.0,
                            skip_checks::Bool=false) where {
                            S <: DiseaseStateSequence,
                            M <: PhyloILM}
    @debug "Initializing $S $M Simulation with the following starting states:" states
    if !skip_checks
      if length(states) != individuals(pop)
        @error "Length of initial disease state vector must match number of individuals"
      elseif !all(in.(states, Ref(convert(DiseaseStates, S))))
        @error "All states in initial disease state vector must be valid within specified epidemic model"
      end
    end
    tr = initialize(TransmissionRates, states, pop, rf, rp)
    rates = initialize(EventRates, tr, states, pop, rf, rp)
    events = Events{S}(states)
    net = TransmissionNetwork(states)
    return new{S, M}(time, 0, pop, rf, rp, sm, copy(states), tr, rates, events, net)
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