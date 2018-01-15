abstract type EventObservations end


"""
SEIR event observations
"""
type SEIR_EventObservations <: EventObservations
  infected::Vector{Float64}
  removed::Vector{Float64}
  individuals::Int64

  function SEIR_EventObservations(infected::Vector{Float64},
                                  removed::Vector{Float64})
    if length(infected) != length(removed)
      error("Infection and removal event time vectors must be of equal length")
    elseif any(infected .>= removed)
      error("All removal observations must occur after an infection observation")
    end
    return new(infected, removed, length(infected))
  end
end


"""
SIR event observations
"""
type SIR_EventObservations <: EventObservations
  infected::Vector{Float64}
  removed::Vector{Float64}
  individuals::Int64

  function SIR_EventObservations(infected::Vector{Float64},
                                 removed::Vector{Float64})
    if length(infected) != length(removed)
      error("Infection and removal event time vectors must be of equal length")
    elseif any(infected .>= removed)
      error("All removal observations must occur after an infection observation")
    end
    return new(infected, removed, length(infected))
  end
end


"""
SEI event observations
"""
type SEI_EventObservations <: EventObservations
  infected::Vector{Float64}
  individuals::Int64

  function SEI_EventObservations(infected::Vector{Float64})
    return new(infected, length(infected))
  end
end


"""
SI event observations
"""
type SI_EventObservations <: EventObservations
  infected::Vector{Float64}
  individuals::Int64

  function SI_EventObservations(infected::Vector{Float64})
    return new(infected, length(infected))
  end
end
