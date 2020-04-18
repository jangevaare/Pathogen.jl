struct Events{T <: EpidemicModel}
  exposure::Union{Nothing, Vector{Float64}}
  infection::Vector{Float64}
  removal::Union{Nothing, Vector{Float64}}
  individuals::Int64

  function Events{T}(e::V, i::V, r::V) where {T <: SEIR, V <: Vector{Float64}}
    if length(unique((length.([e; i; r])))) != 1
      @error "Length of event time vectors must be equal"
    end
    return new{T}(e, i, r, length(i))
  end

  function Events{T}(n::Integer) where T <: SEIR
    return new{T}(fill(NaN, n), fill(NaN, n), fill(NaN, n), n)
  end

  function Events{T}(e::V, i::V) where {T <: SEI, V <: Vector{Float64}}
    if length(unique((length.([e; i])))) != 1
      @error "Length of event time vectors must be equal"
    end
    return new{T}(e, i, nothing, length(i))
  end

  function Events{T}(n::Integer) where T <: SEI
    return new{T}(fill(NaN, n), fill(NaN, n), nothing, n)
  end

  function Events{T}(i::V, r::V) where {T <: SIR, V <: Vector{Float64}}
    if length(unique((length.([i; r])))) != 1
      @error "Length of event time vectors must be equal"
    end
    return new{T}(nothing, i, r, length(i))
  end

  function Events{T}(n::Integer) where T <: SIR
    return new{T}(nothing, fill(NaN, n), fill(NaN, n), n)
  end

  function Events{T}(i::V) where {T <: SI, V <: Vector{Float64}}
    return new{T}(nothing, i, nothing, length(i))
  end

  function Events{T}(n::Integer) where T <: SI
    return new{T}(nothing, fill(NaN, n), nothing, n)
  end
end

function Events{T}(a::Array{Float64,2}) where T <: EpidemicModel
  if size(a, 2) == 3
    return Events{T}(a[:,1], a[:,2], a[:,3])
  elseif size(a, 2) == 2
    return Events{T}(a[:,1], a[:,2])
  elseif size(a, 2) == 1
    return Events{T}(a[:,1])
  else
    @error "Invalid array size for construction of an $(Events{T}) object"
  end
end

function Events{T}(x::Vector{DiseaseState}) where T <: EpidemicModel
  events = Events{T}(length(x))
  for i = 1:length(x)
    if x[i] in _state_progressions[T]
      for j = _state_progressions[T][2:findfirst(Ref(x[i]) .== _state_progressions[T])]
        events[j][i] = -Inf
      end
    else
      @error "Invalid initial state for individual $i"
    end
  end
  return events
end

function Base.show(io::IO, x::Events{T}) where T <: EpidemicModel
  return print(io, "$T model event times (n=$(x.individuals))")
end

function Base.getindex(x::Events{T}, new_state::DiseaseState) where T <: EpidemicModel
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

function Base.getindex(x::Events{T}, states::Vector{DiseaseState}) where T <: EpidemicModel
  y = x[states[1]]
  for i = 2:length(states)
    y = hcat(y, x[states[i]])
  end
  return y
end

function Base.convert(::Type{Array{Float64, 2}}, x::Events{T}) where T <: EpidemicModel
  return x[_state_progressions[T][2:end]]
end

function Base.convert(::Type{Vector{Float64}}, x::Events{T}) where T <: EpidemicModel
  return x[_state_progressions[T][2:end]][:]
end

function Base.convert(::Type{Array{Float64, 2}}, x::Array{Events{T}, 1}) where T <: EpidemicModel
  y = convert(Vector{Float64}, x[1])'
  for i = 2:length(x)
    y = vcat(y, convert(Vector{Float64}, x[i])')
  end
  return y
end

function Base.minimum(x::Events{T}) where T <: EpidemicModel
  y = convert(Array{Float64, 2}, x)
  return minimum(y[.!isnan.(y)])
end

function Base.maximum(x::Events{T}) where T <: EpidemicModel
  y = convert(Array{Float64, 2}, x)
  return maximum(y[.!isnan.(y)])
end

function Statistics.mean(x::Vector{Events{T}}) where T <: EpidemicModel
  return Events{T}(mean([convert(Array{Float64, 2}, i) for i in x]))
end
