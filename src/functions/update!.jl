function update!(events::Events{T},
                 event::Event{T}) where T <: EpidemicModel
  events[event.new_state][event.individual] = event.time
  return events
end

function update!(states::Vector{DiseaseState},
                 event::Event{T}) where T <: EpidemicModel
  states[event.individual] = event.new_state
  return states
end

function update!(net::TransmissionNetwork,
                 xfer::ExogenousTransmission)
  net.external[xfer.individual] = true
  return net
end

function update!(net::TransmissionNetwork,
                 xfer::EndogenousTransmission)
  net.internal[xfer.source, xfer.individual] = true
  return net
end

function update!(net::TransmissionNetwork,
                 xfer::NoTransmission)
  return net
end

function update!(tr::TransmissionRates,
                 event::Event{T},
                 states::Vector{DiseaseState},
                 pop::Population,
                 rf::RiskFunctions{T},
                 rp::RiskParameters{T}) where T <: EpidemicModel
  id = event.individual
  if _new_transmission(event)
    tr.external[id] = 0.0
    tr.internal[:, id] .= 0.0
  end
  if event.new_state == State_I
    transmissibility = rf.transmissibility(rp.transmissibility, pop, id)
    for i in findall(states .== Ref(State_S))
      tr.internal[id, i] = rf.susceptibility(rp.susceptibility, pop, i) *
                           rf.infectivity(rp.infectivity, pop, i, id) *
                           transmissibility
    end
  elseif event.new_state == State_R
    tr.internal[id, :] .= 0.0
  end
  return tr
end

function update!(rates::EventRates{T},
                 tr::TransmissionRates,
                 event::Event{T},
                 states::Vector{DiseaseState},
                 pop::Population,
                 rf::RiskFunctions{T},
                 rp::RiskParameters{T}) where T <: EpidemicModel
  id = event.individual
  if event.new_state == State_E
    rates.exposure[id] = 0.0
    rates.infection[id] = rf.latency(rp.latency, pop, id)
    @logmsg LogLevel(-5000) "Exposure rate total for i = $id updated"  λ = rates.exposure[id]
    @logmsg LogLevel(-5000) "Infection rate for i = $id updated" λ = rates.infection[id]
  elseif event.new_state == State_I
    rates.infection[id] = 0.0
    @logmsg LogLevel(-5000) "Infection rate for i = $id updated" λ = rates.infection[id]
    if T in [SEIR; SIR]
      rates.removal[id] = rf.removal(rp.removal, pop, id)
      @logmsg LogLevel(-5000) "Removal rate for i = $id updated" λ = rates.removal[id]
    end
    if T in [SEIR; SEI]
      @simd for i in findall(states .== Ref(State_S))
        # This assumes `TransmissionRates` already updated!
        rates.exposure[i] += tr.internal[id, i]
        @logmsg LogLevel(-5000) "Exposure rate total for i = $id updated"  λ = rates.exposure[id]
      end
    elseif T in [SIR; SI]
      @simd for i in findall(states .== Ref(State_S))
        # This assumes `TransmissionRates` already updated!
        rates.infection[i] += tr.internal[id, i]
        @logmsg LogLevel(-5000) "Infection rate total for i = $id updated"  λ = rates.infection[id]
      end
    end
  elseif event.new_state == State_R
    rates.removal[id] = 0.0
    @logmsg LogLevel(-5000) "Removal rate for i = $id updated" λ = rates.removal[id]
    if T == SEIR
      @simd for i in findall(states .== Ref(State_S))
        # This assumes `TransmissionRates` already updated!
        rates.exposure[i] = sum(tr.internal[:, i])
        @logmsg LogLevel(-5000) "Exposure rate total for i = $i updated"  λ = rates.exposure[i]
      end
    elseif T == SIR
      @simd for i in findall(states .== Ref(State_S))
        # This assumes `TransmissionRates` already updated!
        rates.infection[i] = sum(tr.internal[:, i])
        @logmsg LogLevel(-5000) "Infection rate for i = $i updated" λ = rates.infection[i]
      end
    end
  end
  return rates
end

function update!(tr::TransmissionRates,
                 rates::EventRates{T},
                 event::Event{T},
                 states::Vector{DiseaseState},
                 pop::Population,
                 rf::RiskFunctions{T},
                 rp::RiskParameters{T}) where T <: EpidemicModel
  update!(tr, event, states, pop, rf, rp)
  update!(rates, tr, event, states, pop, rf, rp)
  return tr, rates
end
