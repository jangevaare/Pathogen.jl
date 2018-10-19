function initialize(::Type{EventRates},
                    tr::TransmissionRates,
                    states::Vector{DiseaseState},
                    pop::Population,
                    rf::RiskFunctions{T},
                    rp::RiskParameters{T}) where T <: EpidemicModel
  n_ids = length(states)
  rates = EventRates{T}(n_ids)
  for i = 1:n_ids
    if states[i] == State_S
      if T in [SEIR; SEI]
        rates.exposure[i] = tr.external[i] + sum(tr.internal[:,i])
      elseif T in [SIR; SI]
        rates.infection[i] = tr.external[i] + sum(tr.internal[:,i])
      end
    elseif states[i] == State_E
      rates.infection[i] = rf.latency(rp.latency, pop, i)
    elseif states[i] == State_I
      if T in [SEIR; SIR]
        rates.removal[i] = rf.removal(rp.removal, pop, i)
      end
    end
  end
  @debug "Initialization of $T EventRates complete" rates = rates[_state_progressions[T][2:end]]
  return rates
end
