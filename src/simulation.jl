"""
simulate!(n::Int64,
          states::States,
          rates::Rates,
          events::Events,
          network::Network,
          population::DataFrame,
          riskfuncs::RiskFunctions,
          riskparams::RiskParameters)

Simulation function
"""
function simulate!(n::Int64,
                   states::States,
                   rates::Rates,
                   events::Events,
                   network::Network,
                   population::DataFrame,
                   riskfuncs::RiskFunctions,
                   riskparams::RiskParameters)
  counter = 0
  time = 0.
  while counter < n && time < Inf
    counter += 1
    time, event = generate_event(rates, time)
    if time < Inf
      update_states!(states, event)
      update_rates!(rates, states, event, population, riskfuncs, riskparams)
      update_events!(events, event, time)
      update_network!(network, event)
    end
  end
  return states, rates, events, network
end


"""
observe(events::Events,
        delay_infected::UnivariateDistribution,
        delay_removed::UnivariateDistribution,
        force = false::Bool)

Make event time observations from a simulation. Force option ensures all
infections are observed.
"""
function observe(events::Events,
                 delay_infected::UnivariateDistribution,
                 delay_removed::UnivariateDistribution,
                 force = false::Bool)
  infected = fill(NaN, events.individuals)
  removed = fill(NaN, events.individuals)
  if force
    @simd for i = 1:events.individuals
      infection_delay_ub = events.removed[i] - events.infected[i]
      infected[i] = events.infected[i] + rand(Truncated(delay_infected, 0., infection_delay_ub))
      removed[i] = events.removed[i] + rand(Truncated(delay_removed, 0., Inf))
    end
  else
    @simd for i = 1:events.individuals
      infection_delay = rand(delay_infected)
      if isnan(events.removed[i]) || infection_delay + events.infected[i] < events.removed[i]
        infected[i] = events.infected[i] + infection_delay
        removed[i] = events.removed[i] + rand(Truncated(delay_removed, 0., Inf))
      end
    end
  end
  return EventObservations(infected, removed)
end


"""
observe(events::Events,
        delay::UnivariateDistribution,
        force = false::Bool)

Make event time observations from a simulation. Force option ensures all
infections are observed.
"""
function observe(events::Events,
                 delay::UnivariateDistribution,
                 force = false::Bool)
  return observe(events, delay, delay, force)
end
