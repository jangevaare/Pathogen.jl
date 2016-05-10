"""
Provides the state of an individual at a specified time
"""
function findstate(events::Events, individual::Int64, time::Float64)
  if !(1 <= individual <= length(events.exposed))
    error("Invalid individual specified")
  end
  if isnan(events.exposed[individual]) || events.exposed[individual] > time
    return "S"
  elseif isnan(events.infected[individual]) || events.infected[individual] > time
    return "E"
  elseif isnan(events.removed[individual]) || events.removed[individual] > time
    return "I"
  elseif events.removed[individual] <= time
    return "R"
  else
    error("Could not determine state of individual")
  end
end


"""
Provides the state of an array of individuals at a specified time
"""
function findstate(events::Events, individuals::Array{Int64}, time::Float64)
  states = fill("", size(individuals))
  for i in eachindex(individuals)
    states[i] = findstate(events, individuals[i], time)
  end
  return states
end


"""
Provides the state of an array of individuals at a specified time
"""
function findstate(events::Events, time::Float64)
  states = fill("", length(events.exposed))
  for i in 1:length(states)
    states[i] = findstate(events, i, time)
  end
  return states
end
