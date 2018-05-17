mutable struct EventObservations{T <: EpidemicModel}
  infected::Vector{Float64}
  removed::Vector{Float64}
  individuals::Int64

  function EventObservations{T}(individuals::Int64) where T <: EpidemicModel
    return _init_EventObservations!(new{T}(), individuals)
  end

  function EventObservations{T}(x::Vector{Float64}, y::Vector{Float64}) where T <: EpidemicModel
    return _init_EventObservations!(new{T}(), x, y)
  end

  function EventObservations{T}(x::Vector{Float64}) where T <: EpidemicModel
    return _init_EventObservations!(new{T}(), x)
  end
end

function _init_EventObservations!(x::EventObservations{T}, individuals::Int64) where T <: Union{SEIR, SIR}
  x.infected = fill(NaN, individuals)
  x.removed = fill(NaN, individuals)
  x.individuals = individuals
  return x
end

function _init_EventObservations!(x::EventObservations{T}, individuals::Int64) where T <: Union{SEI, SI}
  x.infected = fill(NaN, individuals)
  x.individuals = individuals
  return x
end

function _init_EventObservations!(x::EventObservations{T},
                                  infected::Vector{Float64},
                                  removed::Vector{Float64}) where T <: Union{SEIR, SIR}
  if length(infected) != length(removed)
    error("Length of infection and removal times must be equal")
  end
  x.infected = fill(NaN, individuals)
  x.removed = fill(NaN, individuals)
  x.individuals = length(infected)
  return x
end

function _init_EventObservations!(x::EventObservations{T},
                                  infected::Vector{Float64}) where T <: Union{SEI, SI}
  x.infected = fill(NaN, individuals)
  x.individuals = length(infected)
  return x
end
