"""
findstate(events::Events,
          individual::Int64,
          time::Real)

Provides the state of an individual at a specified time
"""
function findstate(events::Events,
                   individual::Int64,
                   time::Real)
  if !(1 <= individual <= events.individuals)
    throw("Invalid individual specified")
  end
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
