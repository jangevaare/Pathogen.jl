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
    for k in findall(states .== Ref(State_I))
      tr.internal[k, i] = rf.susceptibility(rp.susceptibility, pop, i) *
                          rf.transmissibility(rp.transmissibility, pop, i, k) *
                          rf.infectivity(rp.infectivity, pop, k)
    end
  end
  @debug "Initialization of $T TransmissionRates complete" external = tr.external ∑external = sum(tr.external) internal = tr.internal ∑internal = sum(tr.internal)
  return tr
end
