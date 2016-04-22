"""
`Events` type for simulations
"""
type Events
  susceptible::Vector{Float64}
  exposed::Vector{Float64}
  infected::Vector{Float64}
  detected::Vector{Float64}
  removed::Vector{Float64}
  network::Vector{Array{Bool}}
  function Events(population::DataFrame)
    individuals = size(population, 1)
    susceptible = fill(0.0, individuals)
    exposed = fill(NaN, individuals)
    infected = fill(NaN, individuals)
    detected = fill(NaN, individuals)
    removed = fill(NaN, individuals)
    network = Array{Bool}[]
    push!(network, fill(false, individuals))
    push!(network, fill(false, (individuals, individuals)))
    return new(susceptible,
               exposed,
               infected,
               detected,
               removed,
               network)
  end
end

#
# """
# `EventsList` type for tree generation
# """
# type EventsList
#   time::Vector{Float64}
#   event::Vector{Tuple{Int64, Int64}}
#
#   function EventsList(events::Events)
#     return new(time, event)
#   end
#
#   function EventsList(time::Vector{Float64},
#                       event::Vector{Tuple{Int64, Int64}})
#     return new(time, event)
#   end
# end
