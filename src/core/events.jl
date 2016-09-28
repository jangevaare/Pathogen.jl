"""
All epidemic relevant event times
"""
type Events
  exposed::Vector{Float64}
  infected::Vector{Float64}
  removed::Vector{Float64}
  individuals::Int64

  function Events(population::DataFrame)
    individuals = size(population, 1)
    exposed = fill(NaN, individuals)
    infected = fill(NaN, individuals)
    removed = fill(NaN, individuals)
    return new(exposed,
               infected,
               removed,
               individuals)
  end

  function Events(exposed::Vector{Float64},
                  infected::Vector{Float64},
                  removed::Vector{Float64})
    if !(length(exposed) == length(infected) == length(removed))
      throw(BoundsError)
    end
    individuals = length(exposed)
    return new(exposed,
               infected,
               removed,
               individuals)
  end
end
