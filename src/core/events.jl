"""
`Events` type for simulations
"""
type Events
  susceptible::Vector{Float64}
  exposed::Vector{Float64}
  infected::Vector{Float64}
  removed::Vector{Float64}
  network::Vector{Array{Bool}}
  function Events(population::DataFrame)
    individuals = size(population, 1)
    susceptible = fill(0.0, individuals)
    exposed = fill(NaN, individuals)
    infected = fill(NaN, individuals)
    removed = fill(NaN, individuals)
    network = Array{Bool}[]
    push!(network, fill(false, individuals))
    push!(network, fill(false, (individuals, individuals)))
    return new(susceptible,
               exposed,
               infected,
               removed,
               network)
  end
end
