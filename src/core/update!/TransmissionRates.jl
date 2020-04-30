function update!(tr::TransmissionRates,
                 event::Event{S},
                 states::DiseaseStates,
                 pop::Population,
                 rf::RiskFunctions{S},
                 rp::RiskParameters{S};
                 xzero::Bool=false) where S <: DiseaseStateSequence
  id = event.individual
  if _new_transmission(event) && xzero
    tr.external[id] = 0.0
    tr.internal[:, id] .= 0.0
  end
  if event.new_state == State_I
    transmissibility = rf.transmissibility(rp.transmissibility, pop, id)
    @simd for i in findall(states .== Ref(State_S))
      tr.internal[id, i] = rf.susceptibility(rp.susceptibility, pop, i) *
                           rf.infectivity(rp.infectivity, pop, i, id) *
                           transmissibility
    end
  elseif event.new_state == State_R
    tr.internal[id, states .== Ref(State_S)] .= 0.0
  end
  return tr
end

function update!(tr::TransmissionRates,
                 rates::EventRates{S},
                 event::Event{S},
                 states::DiseaseStates,
                 pop::Population,
                 rf::RiskFunctions{S},
                 rp::RiskParameters{S}) where S <: DiseaseStateSequence
  update!(tr, event, states, pop, rf, rp)
  update!(rates, tr, event, states, pop, rf, rp)
  return tr, rates
end