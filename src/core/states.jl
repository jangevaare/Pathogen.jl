"""
State of individuals in a population at a given time
"""
type States
  susceptible::Vector{Bool}
  exposed::Vector{Bool}
  infected::Vector{Bool}
  removed::Vector{Bool}
  individuals::Int64

  function States(population::DataFrame)
    individuals = size(population, 1)
    return new(fill(true, individuals),
               fill(false, individuals),
               fill(false, individuals),
               fill(false, individuals),
               individuals)
  end
end
