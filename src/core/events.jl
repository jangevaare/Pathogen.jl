"""
All epidemic relevant event times
"""
type Events
  exposed::Vector{Float64}
  infected::Vector{Float64}
  removed::Vector{Float64}
  function Events(population::DataFrame)
    individuals = size(population, 1)
    exposed = fill(NaN, individuals)
    infected = fill(NaN, individuals)
    removed = fill(NaN, individuals)
    return new(exposed,
               infected,
               removed)
  end

  function Events(exposed::Vector{Float64},
                  infected::Vector{Float64},
                  removed::Vector{Float64})
    lengths = [length(exposed);
               length(infected);
               length(removed)]
    if !(lengths[1] == lengths[2] == lengths[3])
      error("Event time vectors are not of equal length")
    end
    return new(exposed,
               infected,
               removed)
  end
end


function copy(events::Events)
  return Events(copy(events.exposed),
                copy(events.infected),
                copy(events.removed))
end
