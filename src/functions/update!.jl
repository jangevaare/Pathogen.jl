function update!(events::Events{T},
                 event::Event{T}) where T <: EpidemicModel
  if event.new_state == State_E
    events.exposure[event.individual] = event.time
  elseif event.new_state == State_I
    events.infection[event.individual] = event.time
  elseif event.new_state == State_R
    events.removal[event.individual] = event.time
  end
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
                 pop::DataFrame,
                 rf::RiskFunctions{T},
                 rp::RiskParameters{T}) where T <: EpidemicModel

  id = event.individual
  if _new_transmission(event)
    tr.external[id] = 0.
    tr.internal[:, id] = 0.
    for i in find(states .== State_S)
      tr.internal[id, i] = rf.susceptibility(rp.susceptibility, pop, i) *
                           rf.transmissibility(rp.transmissibility, pop, i, id) *
                           rf.infectivity(rp.infectivity, pop, id)
    end
  elseif event.new_state == State_R
    tr.internal[id, :] = 0.
  end
  return tr
end

function update!(rates::EventRates{T},
                 tr::TransmissionRates,
                 event::Event{T},
                 states::Vector{DiseaseState},
                 pop::DataFrame,
                 rf::RiskFunctions{T},
                 rp::RiskParameters{T}) where T <: EpidemicModel
  id = event.individual
  if event.new_state == State_E
    rates.exposure[id] = 0.
    rates.infection[id] = rf.latency(rp.latency, pop, id)
  elseif event.new_state == State_I
    rates.infection[id] = 0.
    if T in [SEIR; SIR]
      rates.removal[id] = rf.removal(rp.removal, pop, id)
    end
    if T in [SEIR; SEI]
      @simd for i in find(states .== State_S)
        rates.exposure[i] += tr.internal[id, i]
      end
    elseif T in [SIR; SI]
      @simd for i in find(states .== State_S)
        rates.infection[i] += tr.internal[id, i]
      end
    end
  elseif event.new_state == State_R
    rates.removal[id] = 0.
    if T in [SEIR; SEI]
      @simd for i in find(states .== State_S)
        rates.exposure[i] = sum(tr.internal[:, i])
      end
    elseif T in [SIR; SI]
      @simd for i in find(states .== State_S)
        rates.infection[i] = sum(tr.internal[:, i])
      end
    end
  end
  return rates
end
