struct Events{T} where T <: EpidemicModel
  exposed::Vector{Float64}
  infected::Vector{Float64}
  removed::Vector{Float64}
  individuals::Int64

  function Events{T}(individuals::Int64) where T <: EpidemicModel
    return _init_Events!(new(), individuals)
  end

  function Events{T}(individuals::Int64, x::Vector{Float64}, y::Vector{Float64}, z::Vector{Float64}) where T <: EpidemicModel
    return _init_Events!(new(), x, y, z)
  end

  function Events{T}(individuals::Int64, x::Vector{Float64}, y::Vector{Float64}) where T <: EpidemicModel
    return _init_Events!(new(), x, y)
  end

  function Events{T}(individuals::Int64, x::Vector{Float64}) where T <: EpidemicModel
    return _init_Events!(new(), x)
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

function _init_Events!(x::Events{SEIR},
                       exposed::Vector{Float64},
                       infected::Vector{Float64},
                       removed::Vector{Float64})
  if length(unique([length(exposed); length(infected); length(removed)])) !== 1
    error("Event time vectors do not have matching lengths")
  end
  x.exposed = exposed
  x.infected = infected
  x.removed = removed
  x.individuals = length(infected)
  return x
end

function _init_Events!(x::Events{SEI},
                       exposed::Vector{Float64},
                       infected::Vector{Float64})
  if length(unique([length(exposed); length(infected)])) !== 1
    error("Event time vectors do not have matching lengths")
  end
  x.exposed = exposed
  x.infected = infected
  x.individuals = length(infected)
  return x
end

function _init_Events!(x::Events{SIR},
                       infected::Vector{Float64},
                       removed::Vector{Float64})
  if length(unique([length(infected); length(removed)])) !== 1
    error("Event time vectors do not have matching lengths")
  end
  x.infected = infected
  x.removed = removed
  x.individuals = length(infected)
  return x
end

function _init_Events!(x::Events{SI},
                       infected::Vector{Float64})
  x.infected = infected
  x.individuals = length(infected)
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

function copy(events::SI_Events)
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
