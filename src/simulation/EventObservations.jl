struct EventObservations{T <: EpidemicModel}
  infection::Vector{Float64}
  removal::Union{Nothing, Vector{Float64}}
  individuals::Int64

  function EventObservations{T}(i::V, r::V) where {V<:Vector{Float64}, T <: Union{SEIR, SIR}}
    if length(i) != length(r)
      @error "Length of infection and removal times must be equal"
    end
    return new{T}(i, r, length(i))
  end

  function EventObservations{T}(i::V) where {V<:Vector{Float64}, T <: Union{SEI, SI}}
    return new{T}(i, nothing, length(i))
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

function Base.show(io::IO, x::EventObservations{T}) where T <: EpidemicModel
  return print(io, "$T model observations (n=$(x.individuals))")
end
