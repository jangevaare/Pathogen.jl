const DiseaseStates = Vector{DiseaseState}

function DiseaseStates(individuals::Int64)
  return fill(State_S, Int64)
end
