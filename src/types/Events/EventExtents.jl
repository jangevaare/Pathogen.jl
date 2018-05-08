struct EventExtents{T <: EpidemicModel}
  exposure::Float64
  infection::Float64
  removal::Float64

  EventExtents{T}(v...) = _init_EventExtents!(new{T}(), v...)
end

function _init_EventExtents!(x::EventExtents{SEIR}, v...)
  if length(v) !=3
    error("Incorrect number of event extents provided for SEIR models")
  end
  x.exposure = v[1]
  x.infection = v[2]
  x.removal = v[3]
  return x
end

function _init_EventExtents!(x::EventExtents{SEI}, v...)
  if length(v) !=2
    error("Incorrect number of event extents provided for SEI models")
  end
  x.exposure = v[1]
  x.infection = v[2]
  return x
end

function _init_EventExtents!(x::EventExtents{SIR}, v...)
  if length(v) !=2
    error("Incorrect number of event extents provided for SIR models")
  end
  x.infection = v[1]
  x.removal = v[2]
  return x
end

function _init_EventExtents!(x::EventExtents{SI}, v...)
  if length(v) !=1
    error("Incorrect number of event extents provided for SEIR models")
  end
  x.infection = v[1]
  return x
end
