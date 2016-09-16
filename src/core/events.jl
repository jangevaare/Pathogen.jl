"""
All epidemic relevant event times
"""
type Events
  susceptible::Vector{Float64}
  exposed::Vector{Float64}
  infected::Vector{Float64}
  removed::Vector{Float64}
  function Events(population::DataFrame)
    individuals = size(population, 1)
    susceptible = fill(0.0, individuals)
    exposed = fill(NaN, individuals)
    infected = fill(NaN, individuals)
    removed = fill(NaN, individuals)
    return new(susceptible,
               exposed,
               infected,
               removed)
  end
end
