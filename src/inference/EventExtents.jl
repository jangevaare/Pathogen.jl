struct EventExtents{M <: ILM}
  exposure::Union{Nothing, Float64}
  infection::Union{Nothing, Float64}
  removal::Union{Nothing, Float64}

  EventExtents{SEIR}(e, i, r) = new{SEIR}(e, i, r)
  EventExtents{SEI}(e, i)     = new{SEI}(e, i, nothing)
  EventExtents{SIR}(i, r)     = new{SIR}(nothing, i, r)
  EventExtents{SI}(i)         = new{SI}(nothing, i, nothing)
end

function Base.show(io::IO, x::EventExtents{M}) where M <: ILM
  print(io, "$T model event extents")
end
