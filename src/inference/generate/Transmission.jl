function generate(::Type{AbstractTransmission},
                  tr::TransmissionRates,
                  tnd::Nothing,
                  id::Int64)
  external_or_internal = Weights([tr.external[id];
                                  sum(tr.internal[:,id])])
  if sum(external_or_internal) == 0.0
    @error "All transmission rates = 0.0, No transmission can be generated"
    return NoTransmission()
  elseif sample([true; false], external_or_internal)
    @debug "Exogenous tranmission generated"
    return ExogenousTransmission(id)
  else
    source = sample(1:individuals(tr), Weights(tr.internal[:, id]))
    @debug "Endogenous transmission generated (source id = $source)"
    return EndogenousTransmission(id, source)
  end
end

function generate(::Type{AbstractTransmission},
                  tr::TransmissionRates,
                  tnd::TNDistribution,
                  id::Int64)
  external_or_internal = Weights([tr.external[id] * tnd.external[id];
                                 sum(tr.internal[:,id] .* tnd.internal[:,id])])
  if sum(external_or_internal) == 0.0
    @error "All transmission rates = 0.0, No transmission generated for id = $id" external_or_internal
    return NoTransmission()
  elseif sample([true; false], external_or_internal)
    @debug "Exogenous tranmission generated for id = $id" external_or_internal
    return ExogenousTransmission(id)
  else
    source = sample(1:individuals(tr), Weights(tr.internal[:, id] .* tnd.internal[:, id]))
    @debug "Endogenous transmission generated for id = $id (source id = $source)" external_or_internal
    return EndogenousTransmission(id, source)
  end
end