"""
Create and initialize a rate array
"""
function initialize_rates(population::DataFrame,
                          riskfuncs::RiskFunctions,
                          riskparams::RiskParameters)
  individuals = size(population, 1)
  rates = Rates(individuals)
  for i = 1:individuals
    # External exposure
    rates.rates[1][i] = riskfuncs.sparks(riskparams.sparks, population, i)
    # Internal exposure
    for k = 1:individuals
      if i !== k
        rates.rates[2][k, i] = riskfuncs.susceptibility(riskparams.susceptibility, population, i) *
                               riskfuncs.transmissibility(riskparams.transmissibility, population, k) *
                               riskfuncs.infectivity(riskparams.infectivity, population, i, k)
      end
    end
  end
  # Infection onset
  for j = 1:individuals
    rates.rates[3][j] = riskfuncs.latency(riskparams.latency, population, j)
  end
  # Removal
  for k = 1:individuals
    rates.rates[4][k] = riskfuncs.removal(riskparams.removal, population, k)
  end
  # Mask
  rates.mask[1][:] = true
  return rates
end


"""
Initialize a simulation with an events data frame for a phylodynamic individual
level model of infectious disease
"""
function initialize_simulation(population::DataFrame,
                               riskfuncs::RiskFunctions,
                               riskparams::RiskParameters)
  # Initialize rate array
  rates = initialize_rates(population,
                           riskfuncs,
                           riskparams)
  # Initialize events data frame
  events = Events(population)
  # Initialize exposure network
  network = Network(population)
  return rates, events, network
end
