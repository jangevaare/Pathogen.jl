mutable struct EventExtents{T <: EpidemicModel}
  exposure::Float64
  infection::Float64
  removal::Float64

  function EventExtents{T}(e, i, r) where T <: SEIR
    return new{T}(e, i, r)
  end

  function EventExtents{T}(e, i) where T <: SEI
    x = new{T}()
    x.exposure = e
    x.infection = i
    return x
  end

  function EventExtents{T}(i, r) where T <: SIR
    x = new{T}()
    x.infection = i
    x.removal = r
    return x
  end

  function EventExtents{T}(i) where T <: SI
    x = new{T}()
    x.infection = i
    return x
  end
end

function Base.show(io::IO, x::EventExtents{T}) where T <: EpidemicModel
  print(io, "$T model event extents")
end
