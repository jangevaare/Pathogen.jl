struct EventExtents{S <: DiseaseStateSequence}
  exposure::Union{Nothing, Float64}
  infection::Union{Nothing, Float64}
  removal::Union{Nothing, Float64}

  EventExtents{SEIR}(e::F, i::F, r::F) where {F <: Float64} = new{SEIR}(e, i, r)
  EventExtents{SEI}(e::F, i::F) where {F <: Float64}        = new{SEI}(e, i, nothing)
  EventExtents{SIR}(i::F, r::F) where {F <: Float64}        = new{SIR}(nothing, i, r)
  EventExtents{SI}(i::F) where {F <: Float64}               = new{SI}(nothing, i, nothing)
end

function Base.show(io::IO, x::EventExtents{S}) where {S <: DiseaseStateSequence}
  print(io, "$S model event extents")
end
