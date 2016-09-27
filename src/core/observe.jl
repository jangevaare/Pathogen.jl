"""
Event observations
"""
type EventObservations
  infected::Vector{Float64}
  removed::Vector{Float64}
end


"""
Make event time observations from a simulation
"""
function observe(events::Events,
                 delay_infected::UnivariateDistribution,
                 delay_removed::UnivariateDistribution)
  infected = events.infected .+ rand(delay_infected, length(events.infected))
  removed = events.removed .+ rand(delay_removed, length(events.removed))
  return EventObservations(infected, removed)
end


"""
Make event time observations from a simulation
"""
function observe(events::Events,
                 delay::UnivariateDistribution)
  return observe(events, delay, delay)
end
