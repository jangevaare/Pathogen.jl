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
                 pop::DataFrame,
                 rf::RiskFunctions{T},
                 rp::RiskParameters{T};
                 debug_level::Int64=0) where T <: EpidemicModel
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
                 rp::RiskParameters{T};
                 debug_level::Int64=0) where T <: EpidemicModel
  id = event.individual
  if event.new_state == State_E
    rates.exposure[id] = 0.
    rates.infection[id] = rf.latency(rp.latency, pop, id)
    if debug_level >= 5
      println("update!: exposure and infection rates for id = $id are now $(round(rates.exposure[id],3)), and $(round(rates.infection[id],3)) respectively")
    end
  elseif event.new_state == State_I
    rates.infection[id] = 0.
    if debug_level >= 5
      println("update!: infection rate for id = $id is now $(round(rates.infection[id],3))")
    end
    if T in [SEIR; SIR]
      rates.removal[id] = rf.removal(rp.removal, pop, id)
      if debug_level >= 5
        println("update!: removal rate for id = $id is now $(round(rates.removal[id],3))")
      end
    end
    if T in [SEIR; SEI]
      @simd for i in find(states .== State_S)
        # This assumes `TransmissionRates` already updated!
        rates.exposure[i] += tr.internal[id, i]
        if debug_level >= 5
          println("update!: overall exposure rate for id = $i is now $(round(rates.exposure[i], 3))")
        end
      end
    elseif T in [SIR; SI]
      @simd for i in find(states .== State_S)
        # This assumes `TransmissionRates` already updated!
        rates.infection[i] += tr.internal[id, i]
        if debug_level >= 5
          println("update!: overall infection rate for id = $i is now $(round(rates.infection[i], 3))")
        end
      end
    end
  elseif event.new_state == State_R
    rates.removal[id] = 0.
    if debug_level >= 5
      println("update!: removal rate for id = $id is now $(round(rates.removal[id],3))")
    end
    if T == SEIR
      @simd for i in find(states .== State_S)
        # This assumes `TransmissionRates` already updated!
        rates.exposure[i] = sum(tr.internal[:, i])
        if debug_level >= 5
          println("update!: overall exposure rate for id = $i is now $(round(rates.exposure[i], 3))")
        end
      end
    elseif T == SIR
      @simd for i in find(states .== State_S)
        # This assumes `TransmissionRates` already updated!
        rates.infection[i] = sum(tr.internal[:, i])
        if debug_level >= 5
          println("update!: overall infection rate for id = $i is now $(round(rates.infection[i], 3))")
        end
      end
    end
  end
  return rates
end

function update!(tr::TransmissionRates,
                 rates::EventRates{T},
                 event::Event{T},
                 states::Vector{DiseaseState},
                 pop::DataFrame,
                 rf::RiskFunctions{T},
                 rp::RiskParameters{T};
                 debug_level::Int64=0) where T <: EpidemicModel
  update!(tr, event, states, pop, rf, rp)
  update!(rates, tr, event, states, pop, rf, rp)
  return tr, rates
end
