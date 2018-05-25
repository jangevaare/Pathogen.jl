mutable struct Events{T <: EpidemicModel}
  exposure::Vector{Float64}
  infection::Vector{Float64}
  removal::Vector{Float64}
  individuals::Int64

  function Events{T}(n::Int64) where T <: SEIR
    x = new{T}(fill(NaN, n), fill(NaN, n), fill(NaN, n), n)
    return x
  end

  function Events{T}(n::Int64) where T <: SEI
    x = new{T}()
    x.exposure = fill(NaN, n)
    x.infection = fill(NaN, n)
    x.individuals = n
    return x
  end

  function Events{T}(n::Int64) where T <: SIR
    x = new{T}()
    x.infection = fill(NaN, n)
    x.removal = fill(NaN, n)
    x.individuals = n
    return x
  end

  function Events{T}(n::Int64) where T <: SI
    x = new{T}()
    x.infection = fill(NaN, n)
    x.individuals = n
    return x
  end

  function Events{T}(e::V, i::V, r::V) where {T <: SEIR, V <: Vector{Float64}}
    @boundscheck unique(length.(e, i, r)) != 1 && error("Event time vectors of inequal length")
    return new{T}(e, i, r, length(i))
  end

function Events{T}(e::V, i::V) where {T <: SEI, V <: Vector{Float64}}
    @boundscheck unique(length.(e, i)) != 1 && error("Event time vectors of inequal length")
    x = new{T}()
    x.exposure = e
    x.infection = i
    x.individuals = length(i)
    return x
  end

  function Events{T}(i::V, r::V) where {T <: SIR, V <: Vector{Float64}}
    @boundscheck unique(length.(i, r)) != 1 && error("Event time vectors of inequal length")
    x = new{T}()
    x.infection = i
    x.removal = r
    x.individuals = length(i)
    return x
  end

  function Events{T}(i::V) where {T <: SI, V <: Vector{Float64}}
    x = new{T}()
    x.infection = i
    x.individuals = length(i)
    return x
  end
end

function Base.convert(::Type{Array{Float64, 2}}, x::Events{SEIR})
  return [x.exposed x.infected x.removed]
end

function Base.convert(::Type{Array{Float64, 2}}, x::Events{SEI})
  return [x.exposed x.infected]
end

function Base.Base.convert(::Type{Array{Float64, 2}}, x::Events{SIR})
  return [x.infected x.removed]
end

function Base.convert(::Type{Array{Float64, 1}}, x::Events{SI})
  return x.infected
end

function Base.convert(::Type{Array{Float64, 1}}, x::Events{T}) where T <: EpidemicModel
  returnBase.convert(Array{Float64, 2}, x)[:]
end

function Base.convert(::Type{Array{Float64, 2}}, x::Array{Events{T}, 1}) where T <: EpidemicModel
  y = fill(NaN, (length(x), length(convert(Array{Float64, 1}, x[1]))))
  for i = 1:length(x)
    y[i,:] =Base.convert(Array{Float64, 1}, x[i])
  end
  return y
end

function Base.copy(x::Events{SEIR})
  return Events{SEIR}(copy(x.exposed),
                     Base.copy(x.infected),
                     Base.copy(x.removed))
end

function Base.copy(events::Events{SEI})
  return Events{SEI}(copy(x.exposed),
                    Base.copy(x.infected))
end

function Base.copy(events::Events{SIR})
  return Events{SIR}(copy(x.infected),
                    Base.copy(x.removed))
end

function Base.copy(events::Events{SI})
  return Events{SI}(copy(x.infected))
end

function Base.show(io::IO, object::Events{SEIR})
  print(io, "SEIR Events object\nExposures: $(sum(.!isnan.(object.exposed)))\nInfections: $(sum(.!isnan.(object.infected)))\nRemovals: $(sum(.!isnan.(object.removed)))")
end

function Base.show(io::IO, object::Events{SEI})
  print(io, "SEI Events object\nExposures: $(sum(.!isnan.(object.exposed)))\nInfections: $(sum(.!isnan.(object.infected)))")
end

function Base.show(io::IO, object::Events{SIR})
  print(io, "SIR Events object\nInfections: $(sum(.!isnan.(object.infected)))\nRemovals: $(sum(.!isnan.(object.removed)))")
end

function Base.show(io::IO, object::Events{SI})
  print(io, "SI Events object\nInfections: $(sum(.!isnan.(object.infected)))")
end
