function initialize(::Type{TransmissionRates},
                    states::DiseaseStates,
                    pop::Population,
                    rf::RiskFunctions{S},
                    rp::RiskParameters{S}) where {
                    S <: DiseaseStateSequence}
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
  @debug "Initialization of TransmissionRates complete" Σexternal = sum(tr.external) Σinternal = sum(tr.internal)
  return tr
end