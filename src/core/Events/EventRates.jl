struct EventRates{S <: DiseaseStateSequence}
  exposure::Union{Nothing, Vector{Float64}}
  infection::Vector{Float64}
  removal::Union{Nothing, Vector{Float64}}

  function EventRates{S}(n::Int64) where S <: SEIR
    return new{S}(fill(0.0, n), fill(0.0, n), fill(0.0, n))
  end

  function EventRates{S}(n::Int64) where S <: SEI
    return new{S}(fill(0.0, n), fill(0.0, n), nothing)
  end
  
  function EventRates{S}(n::Int64) where S <: SIR
    return new{S}(nothing, fill(0.0, n), fill(0.0, n))
  end
  
  function EventRates{S}(n::Int64) where S <: SI
    return new{S}(nothing, fill(0.0, n), nothing)
  end
end

function individuals(x::EventRates{S}) where {
                     S <: DiseaseStateSequence}
  return length(x.infection)
end

function Base.getindex(x::EventRates{S}, new_state::DiseaseState) where S <: DiseaseStateSequence
  if new_state == State_E && S <: Union{SEIR, SEI}
    return x.exposure
  elseif new_state == State_I
    return x.infection
  elseif new_state == State_R && S <: Union{SEIR, SIR}
    return x.removal
  else
    throw(ErrorException("Invalid indexing disease state"))
  end
end

function Base.getindex(x::EventRates{S}, new_states::DiseaseStates) where S <: DiseaseStateSequence
  y = x[new_states[1]]
  for i=2:length(new_states)
    y = hcat(y, x[new_states[i]])
  end
  return y
end

function Base.sum(x::EventRates{S}) where S <: DiseaseStateSequence
  return sum([sum(x[i]) for i in convert(DiseaseStates, S)[2:end]])
end
