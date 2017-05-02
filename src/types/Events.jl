abstract Events


"""
All SEIR relevant event times
"""
type SEIR_Events <: Events
  exposed::Vector{Float64}
  infected::Vector{Float64}
  removed::Vector{Float64}
  individuals::Int64

  function SEIR_Events(population::DataFrame)
    individuals = size(population, 1)
    exposed = fill(NaN, individuals)
    infected = fill(NaN, individuals)
    removed = fill(NaN, individuals)
    return new(exposed,
               infected,
               removed,
               individuals)
  end

  function SEIR_Events(exposed::Vector{Float64},
                       infected::Vector{Float64},
                       removed::Vector{Float64})
    if !(length(exposed) == length(infected) == length(removed))
      throw(BoundsError)
    end
    individuals = length(infected)
    return new(exposed,
               infected,
               removed,
               individuals)
  end
end


function copy(events::SEIR_Events)
  return SEIR_Events(copy(events.exposed),
                     copy(events.infected),
                     copy(events.removed))
end


"""
All SIR relevant event times
"""
type SIR_Events <: Events
  infected::Vector{Float64}
  removed::Vector{Float64}
  individuals::Int64

  function SIR_Events(population::DataFrame)
    individuals = size(population, 1)
    infected = fill(NaN, individuals)
    removed = fill(NaN, individuals)
    return new(infected,
               removed,
               individuals)
  end

  function SIR_Events(infected::Vector{Float64},
                      removed::Vector{Float64})
    if length(infected) != length(removed)
      throw(BoundsError)
    end
    individuals = length(infected)
    return new(infected,
               removed,
               individuals)
  end
end


function copy(events::SIR_Events)
  return SIR_Events(copy(events.infected),
                    copy(events.removed))
end


"""
All SEI relevant event times
"""
type SEI_Events <: Events
  exposed::Vector{Float64}
  infected::Vector{Float64}
  individuals::Int64

  function SEI_Events(population::DataFrame)
    individuals = size(population, 1)
    exposed = fill(NaN, individuals)
    infected = fill(NaN, individuals)
    return new(exposed,
               infected,
               individuals)
  end

  function SEI_Events(exposed::Vector{Float64},
                      infected::Vector{Float64})
    if length(exposed) != length(infected)
      throw(BoundsError)
    end
    individuals = length(infected)
    return new(exposed,
               infected,
               individuals)
  end
end


function copy(events::SEI_Events)
  return SEI_Events(copy(events.exposed),
                    copy(events.infected))
end


"""
All SI relevant event times
"""
type SI_Events <: Events
  infected::Vector{Float64}
  individuals::Int64

  function SI_Events(population::DataFrame)
    individuals = size(population, 1)
    infected = fill(NaN, individuals)
    return new(infected,
               individuals)
  end

  function SI_Events(infected::Vector{Float64})
    individuals = length(infected)
    return new(infected,
               individuals)
  end
end


function copy(events::SI_Events)
  return SI_Events(copy(events.infected))
end
