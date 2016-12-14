"""
Create and initialize a rate array
"""
function initialize_rates(states::States,
                          population::DataFrame,
                          riskfuncs::RiskFunctions,
                          riskparams::RiskParameters)
  individuals = size(population, 1)
  rates = Rates(individuals)

  # Exposure
  for i in find(states.susceptible)
    # External exposure
    rates.external[i] = riskfuncs.sparks(riskparams.sparks, population, i)

    # Internal exposure
    @simd for k in find(states.infected)
      rates.internal[k, i] = riskfuncs.susceptibility(riskparams.susceptibility, population, i) *
                             riskfuncs.transmissibility(riskparams.transmissibility, population, k) *
                             riskfuncs.infectivity(riskparams.infectivity, population, i, k)
    end
  end

  # Infection onset
  @simd for j in find(states.exposed)
    rates.infection[j] = riskfuncs.latency(riskparams.latency, population, j)
  end

  # Removal
  @simd for k in find(states.infected)
    rates.removal = riskfuncs.removal(riskparams.removal, population, k)
  end

  return rates
end


"""
Initialize a simulation with an events data frame for a phylodynamic individual
level model of infectious disease
"""
function initialize_simulation(population::DataFrame,
                               riskfuncs::RiskFunctions,
                               riskparams::RiskParameters)
  # Initialize state array
  states = States(population)

  # Initialize rate array
  rates = initialize_rates(states,
                           population,
                           riskfuncs,
                           riskparams)

  # Initialize events data frame
  events = Events(population)

  # Initialize exposure network
  network = Network(population)

  return states, rates, events, network
end


"""
A function to generate an event time and event type from a rate array
"""
function generate_event(rates::Rates,
                        time=0.::Float64)
  totals = [sum(rates.external);
            sum(rates.internal);
            sum(rates.infection);
            sum(rates.removal)]
  total = sum(totals)
  if total == Inf
    time = time
    event_type = findfirst(totals .== Inf)
    event_index = findfirst(rates[event_type][:] .== Inf)
  elseif total == 0.
    time = Inf
    event_type = 0
    event_index = 0
  else
    # Generate event time
    time += rand(Exponential(1/total))
    # Generate event type and index
    event_type = findfirst(rand(Multinomial(1, totals/total)))
    event_index = findfirst(rand(Multinomial(1, rates[event_type][:]/totals[event_type])))
  end
  return time, (event_type, event_index)
end


"""
A function to update a `States` object based on an event occurence
"""
function update_states!(states::States,
                        event::Tuple{Int64, Int64})
  # External exposure
  if event[1] == 1
    individual = event[2]
    states.susceptible[individual] = false
    states.exposed[individual] = true

  # Internal exposure
  elseif event[1] == 2
    individual = ind2sub((states.individuals, states.individuals), event[2])[2]
    states.susceptible[individual] = false
    states.exposed[individual] = true

  # Onset of infection
  elseif event[1] == 3
    individual = event[2]
    states.exposed[individual] = false
    states.infected[individual] = true

  # Removal
  elseif event[1] == 4
    individual = event[2]
    states.infected[individual] = false
    states.removed[individual] = true
  end
  return states
end


"""
A function to update a `Rates` object based on an event occurence
"""
function update_rates!(rates::Rates,
                       states::States,
                       event::Tuple{Int64, Int64},
                       population::DataFrame,
                       riskfuncs::RiskFunctions,
                       riskparams::RiskParameters)
  # External exposure
  if event[1] == 1
    individual = event[2]
    rates.external[individual] = 0.
    rates.internal[:, individual] = 0.
    rates.infection[individual] = riskfuncs.latency(riskparams.latency, population, individual)
  # Internal exposure
  elseif event[1] == 2
    individual = ind2sub((rates.individuals, rates.individuals), event[2])[2]
    rates.external[individual] = 0.
    rates.internal[:, individual] = 0.
    rates.infection[individual] = riskfuncs.latency(riskparams.latency, population, individual)
  # Onset of infection
  elseif event[1] == 3
    individual = event[2]
    rates.infection[individual] = 0.
    @simd for i in find(states.susceptible)
      rates.internal[individual, i] = riskfuncs.susceptibility(riskparams.susceptibility, population, i) *
                                      riskfuncs.transmissibility(riskparams.transmissibility, population, individual) *
                                      riskfuncs.infectivity(riskparams.infectivity, population, i, individual)
    end
    rates.removal[individual] = riskfuncs.removal(riskparams.removal, population, individual)
  # Removal
  elseif event[1] == 4
    individual = event[2]
    rates.removal[individual] = 0.
    @simd for i in find(states.susceptible)
      rates.internal[individual, i] = 0.
    end
  end
  return rates
end


"""
A function to update an `Events` object based on an event occurence
"""
function update_events!(events::Events,
                        event::Tuple{Int64, Int64},
                        time::Float64)
  # External exposure
  if event[1] == 1
    individual = event[2]
    events.exposed[individual] = time
  # Internal exposure
  elseif event[1] == 2
    individual = ind2sub((length(events.exposed),length(events.exposed)), event[2])[2]
    events.exposed[individual] = time
  # Onset of infection
  elseif event[1] == 3
    individual = event[2]
    events.infected[individual] = time
  # Removal
  elseif event[1] == 4
    individual = event[2]
    events.removed[individual] = time
  end
  return events
end


"""
A function to update a `Network` object based on an event occurence
"""
function update_network!(network::Network,
                         event::Tuple{Int64, Int64})
  if event[1] == 1
    network.external[event[2]] = true
  elseif event[1] == 2
    network.internal[event[2]] = true
  end
  return network
end


"""
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
Make event time observations from a simulation. Force option ensures all
infections are observed.
"""
function observe(events::Events,
                 delay::UnivariateDistribution,
                 force = false::Bool)
  return observe(events, delay, delay, force)
end
