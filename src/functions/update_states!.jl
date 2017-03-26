"""
update_states!(states::SEIR_States,
               event::Tuple{Int64, Int64})

A function to update a `States` object based on an event occurence
"""
function update_states!(states::SEIR_States,
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
update_states!(states::SIR_States,
               event::Tuple{Int64, Int64})

A function to update a `States` object based on an event occurence
"""
function update_states!(states::SIR_States,
                        event::Tuple{Int64, Int64})
  # External infection
  if event[1] == 1
    individual = event[2]
    states.susceptible[individual] = false
    states.infected[individual] = true

  # Internal infection
  elseif event[1] == 2
    individual = ind2sub((states.individuals, states.individuals), event[2])[2]
    states.susceptible[individual] = false
    states.infected[individual] = true

  # Removal
  elseif event[1] == 3
    individual = event[2]
    states.infected[individual] = false
    states.removed[individual] = true
  end
  return states
end


"""
update_states!(states::SEI_States,
               event::Tuple{Int64, Int64})

A function to update a `States` object based on an event occurence
"""
function update_states!(states::SEI_States,
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
  end
  return states
end


"""
update_states!(states::SI_States,
               event::Tuple{Int64, Int64})

A function to update a `States` object based on an event occurence
"""
function update_states!(states::SI_States,
                        event::Tuple{Int64, Int64})
  # External infection
  if event[1] == 1
    individual = event[2]
    states.susceptible[individual] = false
    states.infected[individual] = true

  # Internal infection
  elseif event[1] == 2
    individual = ind2sub((states.individuals, states.individuals), event[2])[2]
    states.susceptible[individual] = false
    states.infected[individual] = true
  end
  return states
end
