"""
findstate(events::SEIR_Events,
          individual::Int64,
          time::Float64)

Provides the state of an individual at a specified time
"""
function findstate(events::SEIR_Events,
                   individual::Int64,
                   time::Float64)
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


"""
findstate(events::SIR_Events,
          individual::Int64,
          time::Float64)

Provides the state of an individual at a specified time
"""
function findstate(events::SIR_Events,
                   individual::Int64,
                   time::Float64)
  if !(1 <= individual <= events.individuals)
    throw("Invalid individual specified")
  end
  if isnan(events.infected[individual]) || events.infected[individual] > time
    return :S
  elseif isnan(events.removed[individual]) || events.removed[individual] > time
    return :I
  elseif events.removed[individual] <= time
    return :R
  end
end


"""
findstate(events::SEI_Events,
          individual::Int64,
          time::Float64)

Provides the state of an individual at a specified time
"""
function findstate(events::SEI_Events,
                   individual::Int64,
                   time::Float64)
  if !(1 <= individual <= events.individuals)
    throw("Invalid individual specified")
  end
  if isnan(events.exposed[individual]) || events.exposed[individual] > time
    return :S
  elseif isnan(events.infected[individual]) || events.infected[individual] > time
    return :E
  elseif events.infected[individual] <= time
    return :I
  end
end


"""
findstate(events::SI_Events,
          individual::Int64,
          time::Float64)

Provides the state of an individual at a specified time
"""
function findstate(events::SI_Events,
                   individual::Int64,
                   time::Float64)
  if !(1 <= individual <= events.individuals)
    throw("Invalid individual specified")
  end
  if isnan(events.infected[individual]) || events.infected[individual] > time
    return :S
  elseif events.infected[individual] <= time
    return :I
  end
end
