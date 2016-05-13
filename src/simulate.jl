"""
Create an initialize rate array
"""
function initialize_rates(population::DataFrame,
                          risk_funcs::RiskFunctions,
                          risk_params::RiskParameters)
  individuals = size(population, 1)
  rates = Rates(individuals)
  for i = 1:individuals
    # External exposure
    rates.rates[1][i] = risk_funcs.sparks(risk_params.sparks, population, i)
    # Internal exposure
    for k = 1:individuals
      if i !== k
        rates.rates[2][k, i] = risk_funcs.susceptibility(risk_params.susceptibility, population, i) *
                               risk_funcs.transmissibility(risk_params.transmissibility, population, k) *
                               risk_funcs.infectivity(risk_params.infectivity, population, i, k)
      end
    end
  end
  # Infection onset
  for j = 1:individuals
    rates.rates[3][j] = risk_funcs.latency(risk_params.latency, population, j)
  end
  # Detection
  for k = 1:individuals
    rates.rates[4][k] = risk_funcs.detection(risk_params.detection, population, k)
  end
  # Removal
  for k = 1:individuals
    rates.rates[5][k] = risk_funcs.removal(risk_params.removal, population, k)
  end
  # Mask
  rates.mask[1][:] = true
  return rates
end


"""
A function to generate an event time and event type from a rate array
"""
function generate_event(rates::Rates,
                        time=0.::Float64)
  totals = [sum(rates[1]);
            sum(rates[2]);
            sum(rates[3]);
            sum(rates[4]);
            sum(rates[5])]
  total = sum(totals)
  if 0. < total < Inf
    # Generate event time
    time += rand(Exponential(1/total))
    # Generate event type and index
    event_type = findfirst(rand(Multinomial(1, totals/total)))
    event_index = findfirst(rand(Multinomial(1, rates[event_type][:]/totals[event_type])))
  elseif total == Inf
    time = time
    event_type = findfirst(totals .== Inf)
    event_index = findfirst(rates[event_type][:] .== Inf)
  elseif total == 0.
    time = Inf
    event_type = 0
    event_index = 0
  end
  return time, (event_type, event_index)
end


"""
A function to update a `Rates` object based on an event occurence
"""
function update_rates!(rates::Rates,
                       event::Tuple{Int64, Int64})
  # External exposure
  if event[1] == 1
    individual = event[2]
    rates.mask[1][individual] = false
    rates.mask[2][:, individual] = false
    rates.mask[3][individual] = true
  # Internal exposure
  elseif event[1] == 2
    individual = ind2sub(size(rates[2]), event[2])[2]
    rates.mask[1][individual] = false
    rates.mask[2][:, individual] = false
    rates.mask[3][individual] = true
  # Onset of infection
  elseif event[1] == 3
    individual = event[2]
    rates.mask[3][individual] = false
    rates.mask[4][individual] = true
    rates.mask[5][individual] = true
    rates.mask[2][individual, :] = rates.mask[1]
  # Detection
  elseif event[1] == 4
    individual = event[2]
    rates.mask[4][individual] = false
  # Removal
  elseif event[1] == 5
    individual = event[2]
    rates.mask[4][individual] = false
    rates.mask[5][individual] = false
    rates.mask[2][individual, :] = false
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
    events.network[1][event[2]] = true
  # Internal exposure
  elseif event[1] == 2
    individual = ind2sub((length(events.exposed),length(events.exposed)), event[2])[2]
    events.exposed[individual] = time
    events.network[2][event[2]] = true
  # Onset of infection
  elseif event[1] == 3
    individual = event[2]
    events.infected[individual] = time
  # Detection
  elseif event[1] == 4
    individual = event[2]
    events.detected[individual] = time
  # Removal
  elseif event[1] == 5
    individual = event[2]
    events.removed[individual] = time
  end
  return events
end


"""
Initialize a simulation with an events data frame for a phylodynamic individual
level model of infectious disease
"""
function initialize_simulation(population::DataFrame,
                               risk_funcs::RiskFunctions,
                               risk_params::RiskParameters,
                               index_case::Int64)
  # Initialize rate array
  rates = initialize_rates(population,
                           risk_funcs,
                           risk_params)
  # Initialize events data frame
  events = Events(population)
  # Add index case
  update_rates!(rates, (1, index_case))
  update_events!(events, (1, index_case), 0.0)
  return rates, events
end


"""
Initialize a simulation with an events data frame for a phylodynamic individual
level model of infectious disease
"""
function initialize_simulation(population::DataFrame,
                               risk_funcs::RiskFunctions,
                               risk_params::RiskParameters)
  # Initialize rate array
  rates = initialize_rates(population,
                           risk_funcs,
                           risk_params)
  # Initialize events data frame
  events = Events(population)
  return rates, events
end


"""
Simulation function
"""
function simulate!(n::Int64,
                   rates::Rates,
                   events::Events)
  counter = 0
  time = 0.
  while counter < n && time < Inf
    counter += 1
    time, event = generate_event(rates, time)
    if time < Inf
      update_rates!(rates, event)
      update_events!(events, event, time)
    end
  end
  return rates, events
end


"""
Generate phylogenetic tree based on transmission events
"""
function generatetree(events::Events)
  eventtimes = [events.exposed events.detected events.removed]
  eventorder = sortperm(eventtimes[:])
  eventnodes = fill(Nullable{Int64}(), size(eventtimes))
  tree = Tree()
  for i = 1:length(eventorder)
    isnan(eventtimes[eventorder[i]]) && break
    event = ind2sub(size(eventtimes), eventorder[i])
    if event[2] == 1
      # Exposure event
      if events.network[1][event[1]]
        # External exposure
        parentnode = 1
        branch_length = eventtimes[eventorder[i]]
      else
        # Internal exposure
        source = findfirst(events.network[2][:, event[1]])
        priorexposures = events.network[2][source, :][:] & (eventtimes[:, 1] .< eventtimes[eventorder[i]])
        if isnan(eventtimes[source, 2]) || eventtimes[source, 2] > eventtimes[event[1], 1]
          # Undetected exposure source
          if any(priorexposures)
            # Prior exposures from this source
            parentnode = get(eventnodes[priorexposures, 1][indmax(eventtimes[priorexposures, 1])])
            branch_length = eventtimes[eventorder[i]] - maximum(eventtimes[priorexposures, 1])
          else
            # No prior exposures from this source
            parentnode = get(eventnodes[source, 1])
            branch_length = eventtimes[eventorder[i]] - eventtimes[source, 1]
          end
        else
          # Detected exposure source
          if !any(priorexposures) || all(eventtimes[source, 2] .> eventtimes[priorexposures, 1])
            # Detection of exposure source is most recent, relevant event
            parentnode = get(eventnodes[source, 2])
            branch_length = eventtimes[eventorder[i]] - eventtimes[source, 2]
          else
            # Other, prior exposure is most recent, relevant event
            parentnode = get(eventnodes[priorexposures, 1][indmax(eventtimes[priorexposures, 1])])
            branch_length = eventtimes[eventorder[i]] - maximum(eventtimes[priorexposures, 1])
          end
        end
      end
      branch!(tree, parentnode, branch_length)
      eventnodes[eventorder[i]] = length(tree.nodes)
    elseif event[2] == 2
      # Detection event
      priorexposures = events.network[2][event[1], :][:] & (eventtimes[:, 1] .< eventtimes[eventorder[i]])
      if any(priorexposures)
        # Individual has exposed others before detection
        parentnode = get(eventnodes[priorexposures, 1][indmax(eventtimes[priorexposures, 1])])
        branch_length = eventtimes[eventorder[i]] - maximum(eventtimes[priorexposures, 1])
      else
        # Individual has not exposed others before detection
        parentnode = get(eventnodes[event[1], 1])
        branch_length = eventtimes[eventorder[i]] - eventtimes[event[1], 1]
      end
      branch!(tree, parentnode, branch_length)
      eventnodes[eventorder[i]] = length(tree.nodes)
    elseif event[2] == 3
      # Removal event
      priorexposures = events.network[2][event[1], :][:]
      if any(priorexposures)
        # Individual has exposed others prior to removal
        if !isnan(eventtimes[event[1], 2]) && all(eventtimes[priorexposures, 1] .< eventtimes[event[1], 2])
          # Detection is most recent and relevant event
          parentnode = get(eventnodes[event[1], 2])
          branch_length = eventtimes[eventorder[i]] - eventtimes[event[1], 2]
        else
          # Exposure is the most recent and relevant event
          parentnode = get(eventnodes[priorexposures, 1][indmax(eventtimes[priorexposures, 1])])
          branch_length = eventtimes[eventorder[i]] - maximum(eventtimes[priorexposures, 1])
        end
      else
        # No prior exposures
        if !isnan(eventtimes[event[1], 2])
          # Has been previously detected
          parentnode = get(eventnodes[event[1], 2])
          branch_length = eventtimes[eventorder[i]] - eventtimes[event[1], 2]
        else
          # Has not been detected
          parentnode = get(eventnodes[event[1], 1])
          branch_length = eventtimes[eventorder[i]] - eventtimes[event[1], 1]
        end
      end
      branch!(tree, parentnode, branch_length)
      eventnodes[eventorder[i]] = length(tree.nodes)
    end
  end
  return tree, eventnodes[:,2]
end
