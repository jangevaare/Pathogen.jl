function initialize(::Type{EventRates},
                    tr::TransmissionRates,
                    states::DiseaseStates,
                    pop::Population,
                    rf::RiskFunctions{S},
                    rp::RiskParameters{S}) where {
                    S <: DiseaseStateSequence}
  n_ids = length(states)
  rates = EventRates{S}(n_ids)
  for i = 1:n_ids
    if states[i] == State_S
      if S in [SEIR; SEI]
        rates.exposure[i] = tr.external[i] + sum(tr.internal[:,i])
      elseif S in [SIR; SI]
        rates.infection[i] = tr.external[i] + sum(tr.internal[:,i])
      end
    elseif states[i] == State_E
      rates.infection[i] = rf.latency(rp.latency, pop, i)
    elseif states[i] == State_I
      if S in [SEIR; SIR]
        rates.removal[i] = rf.removal(rp.removal, pop, i)
      end
    end
  end
  @debug "Initialization of $T `EventRates` complete" rates = rates[convert(DiseaseStates, S)[2:end]]
  return rates
end