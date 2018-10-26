function sum(x::TransmissionRates)
  return sum(x.external) + sum(x.internal)
end

function sum(x::TransmissionRates, i::Int64)
  if i < 1 | i > x.individuals
    @error "Invalid individual identifier" i
  end
  return x.external[i] + sum(x.internal[:, i])
end

function sum(x::EventRates{T}) where T <: EpidemicModel
  return sum([sum(x[i]) for i in _state_progressions[T][2:end]])
end
