struct EventRates{T <: EpidemicModel}
  exposure::Union{Nothing, Vector{Float64}}
  infection::Union{Nothing, Vector{Float64}}
  removal::Union{Nothing, Vector{Float64}}
  individuals::Int64

  function EventRates{T}(n::Int64) where T <: SEIR
    return new{T}(fill(0.0, n), fill(0.0, n), fill(0.0, n), n)
  end

  function EventRates{T}(n::Int64) where T <: SEI
    return new{T}(fill(0.0, n), fill(0.0, n), nothing, n)
  end
  
  function EventRates{T}(n::Int64) where T <: SIR
    return new{T}(nothing, fill(0.0, n), fill(0.0, n), n)
  end  
  
  function EventRates{T}(n::Int64) where T <: SI
    return new{T}(nothing, fill(0.0, n), nothing, n)
  end
end

function Base.getindex(x::EventRates{T}, new_state::DiseaseState) where T <: EpidemicModel
  if new_state == State_E
    return x.exposure
  elseif new_state == State_I
    return x.infection
  elseif new_state == State_R
    return x.removal
  else
    @error "Unrecognized indexing disease state"
  end
end

function Base.getindex(x::EventRates{T}, new_states::Vector{DiseaseState}) where T <: EpidemicModel
  y = x[new_states[1]]
  for i=2:length(new_states)
    y = hcat(y, x[new_states[i]])
  end
  return y
end

function Base.sum(x::EventRates{T}) where T <: EpidemicModel
  return sum([sum(x[i]) for i in _state_progressions[T][2:end]])
end
