"""
initialize_simulation(population::DataFrame,
                      riskfuncs::SEIR_RiskFunctions,
                      riskparams::SEIR_RiskParameters)

Initialize a simulation with an events data frame for a phylodynamic individual
level model of infectious disease
"""
function initialize_simulation(population::DataFrame,
                               riskfuncs::SEIR_RiskFunctions,
                               riskparams::SEIR_RiskParameters)
  # Initialize state array
  states = SEIR_States(population)

  # Initialize rate array
  rates = initialize_rates(states,
                           population,
                           riskfuncs,
                           riskparams)

  # Initialize events data frame
  events = SEIR_Events(population)

  # Initialize exposure network
  network = Network(population)

  # Generate initial event
  time, event = generate_event(rates)

  # Check to make sure initialization was valid
  if time == Inf
    error("Invalid simulation initialization")
  end

  # Set time to 0.
  time = 0.

  # Update everything
  update_states!(states, event)
  update_rates!(rates, states, event, population, riskfuncs, riskparams)
  update_events!(events, event, time)
  update_network!(network, event)

  return states, rates, events, network
end


"""
initialize_simulation(population::DataFrame,
                      riskfuncs::SIR_RiskFunctions,
                      riskparams::SIR_RiskParameters)

Initialize a simulation with an events data frame for a phylodynamic individual
level model of infectious disease
"""
function initialize_simulation(population::DataFrame,
                               riskfuncs::SIR_RiskFunctions,
                               riskparams::SIR_RiskParameters)
  # Initialize state array
  states = SIR_States(population)

  # Initialize rate array
  rates = initialize_rates(states,
                           population,
                           riskfuncs,
                           riskparams)

  # Initialize events data frame
  events = SIR_Events(population)

  # Initialize exposure network
  network = Network(population)

  # Generate initial event
  time, event = generate_event(rates)

  # Check to make sure initialization was valid
  if time == Inf
    error("Invalid simulation initialization")
  end

  # Set time to 0.
  time = 0.

  # Update everything
  update_states!(states, event)
  update_rates!(rates, states, event, population, riskfuncs, riskparams)
  update_events!(events, event, time)
  update_network!(network, event)

  return states, rates, events, network
end


"""
initialize_simulation(population::DataFrame,
                      riskfuncs::SEI_RiskFunctions,
                      riskparams::SEI_RiskParameters)

Initialize a simulation with an events data frame for a phylodynamic individual
level model of infectious disease
"""
function initialize_simulation(population::DataFrame,
                               riskfuncs::SEI_RiskFunctions,
                               riskparams::SEI_RiskParameters)
  # Initialize state array
  states = SEI_States(population)

  # Initialize rate array
  rates = initialize_rates(states,
                           population,
                           riskfuncs,
                           riskparams)

  # Initialize events data frame
  events = SEI_Events(population)

  # Initialize exposure network
  network = Network(population)

  # Generate initial event
  time, event = generate_event(rates)

  # Check to make sure initialization was valid
  if time == Inf
    error("Invalid simulation initialization")
  end

  # Set time to 0.
  time = 0.

  # Update everything
  update_states!(states, event)
  update_rates!(rates, states, event, population, riskfuncs, riskparams)
  update_events!(events, event, time)
  update_network!(network, event)

  return states, rates, events, network
end


"""
initialize_simulation(population::DataFrame,
                      riskfuncs::SI_RiskFunctions,
                      riskparams::SI_RiskParameters)

Initialize a simulation with an events data frame for a phylodynamic individual
level model of infectious disease
"""
function initialize_simulation(population::DataFrame,
                               riskfuncs::SI_RiskFunctions,
                               riskparams::SI_RiskParameters)
  # Initialize state array
  states = SI_States(population)

  # Initialize rate array
  rates = initialize_rates(states,
                           population,
                           riskfuncs,
                           riskparams)

  # Initialize events data frame
  events = SI_Events(population)

  # Initialize exposure network
  network = Network(population)

  # Generate initial event
  time, event = generate_event(rates)

  # Check to make sure initialization was valid
  if time == Inf
    error("Invalid simulation initialization")
  end

  # Set time to 0.
  time = 0.

  # Update everything
  update_states!(states, event)
  update_rates!(rates, states, event, population, riskfuncs, riskparams)
  update_events!(events, event, time)
  update_network!(network, event)

  return states, rates, events, network
end
