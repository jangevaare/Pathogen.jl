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
