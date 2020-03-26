function update!(rates::EventRates{S},
                 tr::TransmissionRates,
                 event::Event{S},
                 states::DiseaseStates,
                 pop::Population,
                 rf::RiskFunctions{S},
                 rp::RiskParameters{S}) where S <: DiseaseStateSequence
  id = event.individual
  if event.new_state == State_E
    rates.exposure[id] = 0.0
    rates.infection[id] = rf.latency(rp.latency, pop, id)
    @logmsg LogLevel(-5000) "Exposure rate total for i = $id updated"  λ = rates.exposure[id]
    @logmsg LogLevel(-5000) "Infection rate for i = $id updated" λ = rates.infection[id]
  elseif event.new_state == State_I
    rates.infection[id] = 0.0
    @logmsg LogLevel(-5000) "Infection rate for i = $id updated" λ = rates.infection[id]
    if S in [SEIR; SIR]
      rates.removal[id] = rf.removal(rp.removal, pop, id)
      @logmsg LogLevel(-5000) "Removal rate for i = $id updated" λ = rates.removal[id]
    end
    if S in [SEIR; SEI]
      @simd for i in findall(states .== Ref(State_S))
        # This assumes `TransmissionRates` already updated!
        rates.exposure[i] += tr.internal[id, i]
        @logmsg LogLevel(-5000) "Exposure rate total for i = $id updated"  λ = rates.exposure[id]
      end
    elseif S in [SIR; SI]
      @simd for i in findall(states .== Ref(State_S))
        # This assumes `TransmissionRates` already updated!
        rates.infection[i] += tr.internal[id, i]
        @logmsg LogLevel(-5000) "Infection rate total for i = $id updated"  λ = rates.infection[id]
      end
    end
  elseif event.new_state == State_R
    rates.removal[id] = 0.0
    @logmsg LogLevel(-5000) "Removal rate for i = $id updated" λ = rates.removal[id]
    if S ==SEIR
      @simd for i in findall(states .== Ref(State_S))
        # This assumes `TransmissionRates` already updated!
        rates.exposure[i] = sum(tr.internal[:, i])
        @logmsg LogLevel(-5000) "Exposure rate total for i = $i updated"  λ = rates.exposure[i]
      end
    elseif S ==SIR
      @simd for i in findall(states .== Ref(State_S))
        # This assumes `TransmissionRates` already updated!
        rates.infection[i] = sum(tr.internal[:, i])
        @logmsg LogLevel(-5000) "Infection rate for i = $i updated" λ = rates.infection[i]
      end
    end
  end
  return rates
end