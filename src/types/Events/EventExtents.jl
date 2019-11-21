struct EventExtents{T <: EpidemicModel}
  exposure::Union{Nothing, Real}
  infection::Union{Nothing, Real}
  removal::Union{Nothing, Real}

  EventExtents{SEIR}(e, i, r) = new{SEIR}(e, i, r)
  EventExtents{SEI}(e, i)     = new{SEI}(e, i, nothing)
  EventExtents{SIR}(i, r)     = new{SIR}(nothing, i, r)
  EventExtents{SI}(i)         = new{SI}(nothing, i, nothing)
end

function Base.show(io::IO, x::EventExtents{T}) where T <: EpidemicModel
  print(io, "$T model event extents")
end
