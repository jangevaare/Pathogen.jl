struct EventObservations{T <: DiseaseStateSequence}
  infection::Vector{Float64}
  removal::Union{Nothing, Vector{Float64}}

  function EventObservations{T}(i::V, r::V) where {V<:Vector{Float64}, T <: Union{SEIR, SIR}}
    if length(i) != length(r)
      error("Length of infection and removal times must be equal")
    end
    return new{T}(i, r)
  end

  function EventObservations{T}(i::V) where {V<:Vector{Float64}, T <: Union{SEI, SI}}
    return new{T}(i, nothing)
  end
end

function EventObservations{T}(i::Array{Float64, 2}) where T <: Union{SEI, SI}
  if size(i, 2) != 1
    @error "Invalid Array dimensions for observations of a $T model"
  end
  return EventObservations{T}(i[:,1])
end

function EventObservations{T}(ir::Array{Float64, 2}) where T <: Union{SEIR, SIR}
  if size(ir, 2) != 2
    @error "Invalid Array dimensions for observations of a $T model"
  end
  return EventObservations(ir[:, 1], ir[:, 2])
end

function individuals(x::EventObservations{M}) where {
                     M <: DiseaseStateSequence}
  return length(x.infection)
end

function Base.show(io::IO, x::EventObservations{T}) where T <: DiseaseStateSequence
  return print(io, "$T model observations (n=$(individuals(x)))")
end

function Base.getindex(x::EventObservations{T}, new_state::DiseaseState) where T <: DiseaseStateSequence
  if new_state == State_I
    return x.infection
  elseif (new_state == State_R) && (T ∈ [SEIR, SIR])
    return x.removal
  else
    error("Unrecognized indexing disease state")
  end
end

function Base.getindex(x::EventObservations{T}, states::DiseaseStates) where T <: DiseaseStateSequence
  y = x[states[1]]
  for i = 2:length(states)
    y = hcat(y, x[states[i]])
  end
  return y
end

function Base.convert(::Type{Array{Float64, 2}}, x::EventObservations{T}) where T <: DiseaseStateSequence
  states = T ∈ [SEI, SI] ? [State_I] : [State_I, State_R]
  return x[states]
end

function Base.convert(::Type{Vector{Float64}}, x::EventObservations{T}) where T <: DiseaseStateSequence
  states = T ∈ [SEI, SI] ? [State_I] : [State_I, State_R]
  return x[states][:]
end

function Base.minimum(x::EventObservations{T}) where T <: DiseaseStateSequence
  y = convert(Array{Float64, 2}, x)
  return minimum(y[.!isnan.(y)])
end

function Base.maximum(x::EventObservations{T}) where T <: DiseaseStateSequence
  y = convert(Array{Float64, 2}, x)
  return maximum(y[.!isnan.(y)])
end
