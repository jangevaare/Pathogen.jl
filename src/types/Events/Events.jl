mutable struct Events{T <: EpidemicModel}
  exposed::Vector{Float64}
  infected::Vector{Float64}
  removed::Vector{Float64}
  individuals::Int64

  function Events{T}(individuals::Int64) where T <: EpidemicModel
    return _init_Events!(new{T}(), individuals)
  end

  function Events{T}(v...) where T <: EpidemicModel
    function _init_Events!(x::Events{SIR}, v)
      if unique(length.(v)) != 1
        error("Mismatch in length of event time vectors")
      end
      x.infected = v[1]
      x.removed = v[2]
      x.individuals = length(x.infected)
      return x
    end
    return _init_Events!(new{T}(), v)
  end
end

function _init_Events!(x::Events{SEIR}, individuals::Int64)
  x.exposed = fill(NaN, individuals)
  x.infected = fill(NaN, individuals)
  x.removed = fill(NaN, individuals)
  x.individuals = individuals
  return x
end

function _init_Events!(x::Events{SEI}, individuals::Int64)
  x.exposed = fill(NaN, individuals)
  x.infected = fill(NaN, individuals)
  x.individuals = individuals
  return x
end

function _init_Events!(x::Events{SIR}, individuals::Int64)
  x.infected = fill(NaN, individuals)
  x.removed = fill(NaN, individuals)
  x.individuals = individuals
  return x
end

function _init_Events!(x::Events{SI}, individuals::Int64)
  x.infected = fill(NaN, individuals)
  x.individuals = individuals
  return x
end

function _init_Events!(x::Events{SEIR}, v)
  if length(v) !=3
    error("Incorrect number of event time vectors provided for SEIR models")
  end
  x.exposed = v[1]
  x.infected = v[2]
  x.removed = v[3]
  x.individuals = length(x.infected)
  return x
end

function _init_Events!(x::Events{SEI}, v)
  if length(v) !=2
    error("Incorrect number of event time vectors provided for SEI models")
  end
  x.exposed = v[1]
  x.infected = v[2]
  x.individuals = length(x.infected)
  return x
end

function _init_Events!(x::Events{SIR}, v)
  if length(v) !=2
    error("Incorrect number of event time vectors provided for SIR models")
  end
  x.infected = v[1]
  x.removed = v[2]
  x.individuals = length(x.infected)
  return x
end

function _init_Events!(x::Events{SI}, v)
  if length(v) !=1
    error("Incorrect number of event time vectors provided for SI models")
  end
  x.infected = v[1]
  x.individuals = length(x.infected)
  return x
end

function convert(::Type{Array{Float64, 2}}, x::Events{SEIR})
  return [x.exposed x.infected x.removed]
end

function convert(::Type{Array{Float64, 2}}, x::Events{SEI})
  return [x.exposed x.infected]
end

function convert(::Type{Array{Float64, 2}}, x::Events{SIR})
  return [x.infected x.removed]
end

function convert(::Type{Array{Float64, 1}}, x::Events{SI})
  return x.infected
end

function convert(::Type{Array{Float64, 1}}, x::Events{T}) where T <: EpidemicModel
  return convert(Array{Float64, 2}, x)[:]
end

function convert(::Type{Array{Float64, 2}}, x::Array{Events{T}, 1}) where T <: EpidemicModel
  y = fill(NaN, (length(x), length(convert(Array{Float64, 1}, x[1]))))
  for i = 1:length(x)
    y[i,:] = convert(Array{Float64, 1}, x[i])
  end
  return y
end

function copy(x::Events{SEIR})
  return Events{SEIR}(copy(x.exposed),
                      copy(x.infected),
                      copy(x.removed))
end

function copy(events::Events{SEI})
  return Events{SEI}(copy(x.exposed),
                     copy(x.infected))
end

function copy(events::Events{SIR})
  return Events{SIR}(copy(x.infected),
                     copy(x.removed))
end

function copy(events::Events{SI})
  return Events{SI}(copy(x.infected))
end

function show(io::IO, object::Events{SEIR})
  print(io, "SEIR Events object\nExposures: $(sum(.!isnan.(object.exposed)))\nInfections: $(sum(.!isnan.(object.infected)))\nRemovals: $(sum(.!isnan.(object.removed)))")
end

function show(io::IO, object::Events{SEI})
  print(io, "SEI Events object\nExposures: $(sum(.!isnan.(object.exposed)))\nInfections: $(sum(.!isnan.(object.infected)))")
end

function show(io::IO, object::Events{SIR})
  print(io, "SIR Events object\nInfections: $(sum(.!isnan.(object.infected)))\nRemovals: $(sum(.!isnan.(object.removed)))")
end

function show(io::IO, object::Events{SI})
  print(io, "SI Events object\nInfections: $(sum(.!isnan.(object.infected)))")
end
