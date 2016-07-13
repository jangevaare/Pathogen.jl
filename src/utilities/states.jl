"""
Provides the state of an individual at a specified time
"""
function findstate(events::Events, individual::Int64, time::Float64, number=false::Bool)
  if !(1 <= individual <= length(events.susceptible))
    error("Invalid individual specified")
  end
  if number
    if isnan(events.exposed[individual]) || events.exposed[individual] > time
      return 1
    elseif isnan(events.infected[individual]) || events.infected[individual] > time
      return 2
    elseif isnan(events.removed[individual]) || events.removed[individual] > time
      return 3
    elseif events.removed[individual] <= time
      return 4
    end
  else
    if isnan(events.exposed[individual]) || events.exposed[individual] > time
      return :S
    elseif isnan(events.infected[individual]) || events.infected[individual] > time
      return :E
    elseif isnan(events.removed[individual]) || events.removed[individual] > time
      return :I
    elseif events.removed[individual] <= time
      return :R
    end
  end
end


"""
Provides the state of an array of individuals at a specified time
"""
function findstate(events::Events, individuals::Array{Int64}, time::Float64, number=false::Bool)
  if number
    states = fill(0, size(individuals))
  else
    states = fill("", size(individuals))
  end
  for i in eachindex(individuals)
    states[i] = findstate(events, individuals[i], time, number)
  end
  return states
end


"""
Provides the state of all individuals at a specified time
"""
function findstate(events::Events, time::Float64, number=false::Bool)
  if number
    states = Int64[]
  else
    states = Symbol[]
  end
  for i in 1:length(events.susceptible)
    push!(states, findstate(events, i, time, number))
  end
  return states
end
