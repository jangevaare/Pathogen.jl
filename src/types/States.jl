abstract type States end


"""
State of individuals in a population at a given time in an SEIR model
"""
type SEIR_States <: States
  susceptible::Vector{Bool}
  exposed::Vector{Bool}
  infected::Vector{Bool}
  removed::Vector{Bool}
  individuals::Int64

  function SEIR_States(population::DataFrame)
    individuals = size(population, 1)
    return new(fill(true, individuals),
               fill(false, individuals),
               fill(false, individuals),
               fill(false, individuals),
               individuals)
  end
end


"""
State of individuals in a population at a given time in an SIR model
"""
type SIR_States <: States
  susceptible::Vector{Bool}
  infected::Vector{Bool}
  removed::Vector{Bool}
  individuals::Int64

  function SIR_States(population::DataFrame)
    individuals = size(population, 1)
    return new(fill(true, individuals),
               fill(false, individuals),
               fill(false, individuals),
               individuals)
  end
end


"""
State of individuals in a population at a given time in an SEI model
"""
type SEI_States <: States
  susceptible::Vector{Bool}
  exposed::Vector{Bool}
  infected::Vector{Bool}
  individuals::Int64

  function SEI_States(population::DataFrame)
    individuals = size(population, 1)
    return new(fill(true, individuals),
               fill(false, individuals),
               fill(false, individuals),
               individuals)
  end
end


"""
State of individuals in a population at a given time in an SI model
"""
type SI_States <: States
  susceptible::Vector{Bool}
  infected::Vector{Bool}
  individuals::Int64

  function SI_States(population::DataFrame)
    individuals = size(population, 1)
    return new(fill(true, individuals),
               fill(false, individuals),
               individuals)
  end
end
