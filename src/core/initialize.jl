function initialize(::Type{TransmissionRates},
                    states::Vector{DiseaseState},
                    pop::Population,
                    rf::RiskFunctions{T},
                    rp::RiskParameters{T}) where T <: EpidemicModel
  n_ids = length(states)
  tr = TransmissionRates(n_ids)
  for i in findall(states .== Ref(State_S))
    # External exposure
    #tr.external[i] = rf.susceptibility(rp.susceptibility, pop, i) * rf.sparks(rp.sparks, pop, i)
    tr.external[i] = rf.sparks(rp.sparks, pop, i)
    # Internal exposure
    susceptibility = rf.susceptibility(rp.susceptibility, pop, i)
    for k in findall(states .== Ref(State_I))
      tr.internal[k, i] = susceptibility *
                          rf.infectivity(rp.infectivity, pop, i, k) *
                          rf.transmissibility(rp.transmissibility, pop, k)
    end
  end
  @debug "Initialization of $T TransmissionRates complete" ∑external = sum(tr.external) ∑internal = sum(tr.internal)
  return tr
end

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
        rates.exposure[i] = sum(tr, i)
      elseif T in [SIR; SI]
        rates.infection[i] = sum(tr, i)
      end
    elseif states[i] == State_E
      rates.infection[i] = rf.latency(rp.latency, pop, i)
    elseif states[i] == State_I
      if T in [SEIR; SIR]
        rates.removal[i] = rf.removal(rp.removal, pop, i)
      end
    end
  end
  @debug "Initialization of $T EventRates complete" rates = [sum(rates[k]) for k in _state_progressions[T][2:end]]
  return rates
end