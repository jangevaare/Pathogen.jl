function initialize(::Type{TransmissionRates},
                    states::Vector{DiseaseState},
                    pop::DataFrame,
                    rf::RiskFunctions{T},
                    rp::RiskParameters{T}) where T <: EpidemicModel

  n_ids = length(states)
  tr = TransmissionRates(n_ids)

  for i in find(states .== State_S)
    # External exposure
    #tr.external[i] = rf.susceptibility(rp.susceptibility, pop, i) * rf.sparks(rp.sparks, pop, i)
    tr.external[i] = rf.sparks(rp.sparks, pop, i)

    # Internal exposure
    @simd for k in find(states .== State_I)
      tr.internal[k, i] = rf.susceptibility(rp.susceptibility, pop, i) *
                          rf.transmissibility(rp.transmissibility, pop, k) *
                          rf.infectivity(rp.infectivity, pop, i, k)
    end
  end
  return tr
end

function initialize(::Type{EventRates{T}},
                    tr::TransmissionRates,
                    states::Vector{DiseaseState},
                    pop::DataFrame,
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
  return rates
end

function initialize(::Type{Trace}, mcmc::MCMC{T}) where T <: EpidemicModel
  events = generate(Events, mcmc)
  rparams = generate(RiskParameters, mcmc)
  #TODO
end
