"""
update_events!(events::SEIR_Events,
               event::SEIR_Event,
               time::Float64)

A function to update an `Events` object based on an event occurence
"""
function update_events!(events::SEIR_Events,
                        event::SEIR_Event,
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
update_events!(events::SIR_Events,
               event::SIR_Event,
               time::Float64)

A function to update an `Events` object based on an event occurence
"""
function update_events!(events::SIR_Events,
                        event::SIR_Event,
                        time::Float64)
  # External infection
  if event[1] == 1
    individual = event[2]
    events.infected[individual] = time
  # Internal infection
  elseif event[1] == 2
    individual = ind2sub((length(events.infected),length(events.infected)), event[2])[2]
    events.infected[individual] = time
  # Removal
  elseif event[1] == 3
    individual = event[2]
    events.removed[individual] = time
  end
  return events
end


"""
update_events!(events::SEI_Events,
               event::SEI_Event,
               time::Float64)

A function to update an `Events` object based on an event occurence
"""
function update_events!(events::SEI_Events,
                        event::SEI_Event,
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
  end
  return events
end


"""
update_events!(events::SI_Events,
               event::SI_Event,
               time::Float64)

A function to update an `Events` object based on an event occurence
"""
function update_events!(events::SI_Events,
                        event::SI_Event,
                        time::Float64)
  # External infection
  if event[1] == 1
    individual = event[2]
    events.infected[individual] = time
  # Internal infection
  elseif event[1] == 2
    individual = ind2sub((length(events.infected),length(events.infected)), event[2])[2]
    events.infected[individual] = time
  end
  return events
end
