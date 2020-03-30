struct Events{S <: DiseaseStateSequence}
  exposure::Union{Nothing, Vector{Float64}}
  infection::Vector{Float64}
  removal::Union{Nothing, Vector{Float64}}
  individuals::Int64

  function Events{S}(e::V, i::V, r::V) where {S <: SEIR, V <: Vector{Float64}}
    if length(unique((length.([e; i; r])))) != 1
      @error "Length of event time vectors must be equal"
    end
    return new{S}(e, i, r, length(i))
  end

  function Events{S}(n::Int64) where S <: SEIR
    return new{S}(fill(NaN, n), fill(NaN, n), fill(NaN, n), n)
  end

  function Events{S}(e::V, i::V) where {S <: SEI, V <: Vector{Float64}}
    if length(unique((length.([e; i])))) != 1
      @error "Length of event time vectors must be equal"
    end
    return new{S}(e, i, nothing, length(i))
  end

  function Events{S}(n::Int64) where S <: SEI
    return new{S}(fill(NaN, n), fill(NaN, n), nothing, n)
  end

  function Events{S}(i::V, r::V) where {S <: SIR, V <: Vector{Float64}}
    if length(unique((length.([i; r])))) != 1
      @error "Length of event time vectors must be equal"
    end
    return new{S}(nothing, i, r, length(i))
  end

  function Events{S}(n::Int64) where S <: SIR
    return new{S}(nothing, fill(NaN, n), fill(NaN, n), n)
  end

  function Events{S}(i::V) where {S <: SI, V <: Vector{Float64}}
    return new{S}(nothing, i, nothing, length(i))
  end

  function Events{S}(n::Int64) where S <: SI
    return new{S}(nothing, fill(NaN, n), nothing, n)
  end
end

function Events{S}(a::Array{Float64,2}) where S <: DiseaseStateSequence
  if size(a, 2) == 3
    return Events{S}(a[:,1], a[:,2], a[:,3])
  elseif size(a, 2) == 2
    return Events{S}(a[:,1], a[:,2])
  elseif size(a, 2) == 1
    return Events{S}(a[:,1])
  else
    @error "Invalid array size for construction of an $(Events{S}) object"
  end
end

function Events{S}(x::DiseaseStates) where S <: DiseaseStateSequence
  events = Events{S}(length(x))
  for i = 1:length(x)
    if x[i] in convert(DiseaseStates, S)
      for j = convert(DiseaseStates, S)[2:findfirst(Ref(x[i]) .== convert(DiseaseStates, S))]
        events[j][i] = -Inf
      end
    else
      @error "Invalid initial state for individual $i"
    end
  end
  return events
end

function Base.show(io::IO, x::Events{S}) where S <: DiseaseStateSequence
  return print(io, "$S model event times (n=$(x.individuals))")
end

function Base.getindex(x::Events{S}, new_state::DiseaseState) where S <: DiseaseStateSequence
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

function Base.getindex(x::Events{S}, states::DiseaseStates) where S <: DiseaseStateSequence
  y = x[states[1]]
  for i = 2:length(states)
    y = hcat(y, x[states[i]])
  end
  return y
end

function Base.convert(::Type{Array{Float64, 2}}, x::Events{S}) where S <: DiseaseStateSequence
  return x[convert(DiseaseStates, S)[2:end]]
end

function Base.convert(::Type{Vector{Float64}}, x::Events{S}) where S <: DiseaseStateSequence
  return x[convert(DiseaseStates, S)[2:end]][:]
end

function Base.convert(::Type{Array{Float64, 2}}, x::Array{Events{S}, 1}) where S <: DiseaseStateSequence
  y = convert(Vector{Float64}, x[1])'
  for i = 2:length(x)
    y = vcat(y, convert(Vector{Float64}, x[i])')
  end
  return y
end

function Base.minimum(x::Events{S}) where S <: DiseaseStateSequence
  y = convert(Array{Float64, 2}, x)
  return minimum(y[.!isnan.(y)])
end

function Base.maximum(x::Events{S}) where S <: DiseaseStateSequence
  y = convert(Array{Float64, 2}, x)
  return maximum(y[.!isnan.(y)])
end

function Statistics.mean(x::Vector{Events{S}}) where S <: DiseaseStateSequence
  return Events{S}(mean([convert(Array{Float64, 2}, i) for i in x]))
end
