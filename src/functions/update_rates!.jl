"""
update_rates!(rates::SEIR_Rates,
              states::SEIR_States,
              event::SEIR_Event,
              population::DataFrame,
              riskfuncs::SEIR_RiskFunctions,
              riskparams::SEIR_RiskParameters)

A function to update a `Rates` object based on an event occurence
"""
function update_rates!(rates::SEIR_Rates,
                       states::SEIR_States,
                       event::SEIR_Event,
                       population::DataFrame,
                       riskfuncs::SEIR_RiskFunctions,
                       riskparams::SEIR_RiskParameters)
  # External exposure
  if event[1] == 1
    individual = event[2]
    rates.exposure.external[individual] = 0.
    rates.exposure.internal[:, individual] = 0.
    rates.infection[individual] = riskfuncs.latency(riskparams.latency, population, individual)
  # Internal exposure
  elseif event[1] == 2
    individual = ind2sub((rates.individuals, rates.individuals), event[2])[2]
    rates.exposure.external[individual] = 0.
    rates.exposure.internal[:, individual] = 0.
    rates.infection[individual] = riskfuncs.latency(riskparams.latency, population, individual)
  # Onset of infection
  elseif event[1] == 3
    individual = event[2]
    rates.infection[individual] = 0.
    @simd for i in find(states.susceptible)
      rates.exposure.internal[individual, i] = riskfuncs.susceptibility(riskparams.susceptibility, population, i) *
                                               riskfuncs.transmissibility(riskparams.transmissibility, population, individual) *
                                               riskfuncs.infectivity(riskparams.infectivity, population, i, individual)
    end
    rates.removal[individual] = riskfuncs.removal(riskparams.removal, population, individual)
  # Removal
  elseif event[1] == 4
    individual = event[2]
    rates.removal[individual] = 0.
    @simd for i in find(states.susceptible)
      rates.exposure.internal[individual, i] = 0.
    end
  end
  return rates
end


"""
update_rates!(rates::SIR_Rates,
              states::SIR_States,
              event::SIR_Event,
              population::DataFrame,
              riskfuncs::SIR_RiskFunctions,
              riskparams::SIR_RiskParameters)

A function to update a `Rates` object based on an event occurence
"""
function update_rates!(rates::SIR_Rates,
                       states::SIR_States,
                       event::SIR_Event,
                       population::DataFrame,
                       riskfuncs::SIR_RiskFunctions,
                       riskparams::SIR_RiskParameters)
  # External infection
  if event[1] == 1
    individual = event[2]
    rates.infection.external[individual] = 0.
    rates.infection.internal[:, individual] = 0.
    @simd for i in find(states.susceptible)
      rates.infection.internal[individual, i] = riskfuncs.susceptibility(riskparams.susceptibility, population, i) *
                                                riskfuncs.transmissibility(riskparams.transmissibility, population, individual) *
                                                riskfuncs.infectivity(riskparams.infectivity, population, i, individual)
    end
    rates.removal[individual] = riskfuncs.removal(riskparams.removal, population, individual)
  # Internal infection
  elseif event[1] == 2
    individual = ind2sub((rates.individuals, rates.individuals), event[2])[2]
    rates.infection.external[individual] = 0.
    rates.infection.internal[:, individual] = 0.
    @simd for i in find(states.susceptible)
      rates.infection.internal[individual, i] = riskfuncs.susceptibility(riskparams.susceptibility, population, i) *
                                                riskfuncs.transmissibility(riskparams.transmissibility, population, individual) *
                                                riskfuncs.infectivity(riskparams.infectivity, population, i, individual)
    end
    rates.removal[individual] = riskfuncs.removal(riskparams.removal, population, individual)
  # Removal
  elseif event[1] == 3
    individual = event[2]
    rates.removal[individual] = 0.
    @simd for i in find(states.susceptible)
      rates.infection.internal[individual, i] = 0.
    end
  end
  return rates
end


"""
update_rates!(rates::SEI_Rates,
              states::SEI_States,
              event::SEI_Event,
              population::DataFrame,
              riskfuncs::SEI_RiskFunctions,
              riskparams::SEI_RiskParameters)

A function to update a `Rates` object based on an event occurence
"""
function update_rates!(rates::SEI_Rates,
                       states::SEI_States,
                       event::SEI_Event,
                       population::DataFrame,
                       riskfuncs::SEI_RiskFunctions,
                       riskparams::SEI_RiskParameters)
  # External exposure
  if event[1] == 1
    individual = event[2]
    rates.exposure.external[individual] = 0.
    rates.exposure.internal[:, individual] = 0.
    rates.infection[individual] = riskfuncs.latency(riskparams.latency, population, individual)
  # Internal exposure
  elseif event[1] == 2
    individual = ind2sub((rates.individuals, rates.individuals), event[2])[2]
    rates.exposure.external[individual] = 0.
    rates.exposure.internal[:, individual] = 0.
    rates.infection[individual] = riskfuncs.latency(riskparams.latency, population, individual)
  # Onset of infection
  elseif event[1] == 3
    individual = event[2]
    rates.infection[individual] = 0.
    @simd for i in find(states.susceptible)
      rates.exposure.internal[individual, i] = riskfuncs.susceptibility(riskparams.susceptibility, population, i) *
                                               riskfuncs.transmissibility(riskparams.transmissibility, population, individual) *
                                               riskfuncs.infectivity(riskparams.infectivity, population, i, individual)
    end
  end
  return rates
end


"""
update_rates!(rates::SI_Rates,
              states::SI_States,
              event::SI_Event,
              population::DataFrame,
              riskfuncs::SI_RiskFunctions,
              riskparams::SI_RiskParameters)

A function to update a `Rates` object based on an event occurence
"""
function update_rates!(rates::SI_Rates,
                       states::SI_States,
                       event::SI_Event,
                       population::DataFrame,
                       riskfuncs::SI_RiskFunctions,
                       riskparams::SI_RiskParameters)
  # External infection
  if event[1] == 1
    individual = event[2]
    rates.infection.external[individual] = 0.
    rates.infection.internal[:, individual] = 0.
    @simd for i in find(states.susceptible)
      rates.infection.internal[individual, i] = riskfuncs.susceptibility(riskparams.susceptibility, population, i) *
                                                riskfuncs.transmissibility(riskparams.transmissibility, population, individual) *
                                                riskfuncs.infectivity(riskparams.infectivity, population, i, individual)
    end
  # Internal infection
  elseif event[1] == 2
    individual = ind2sub((rates.individuals, rates.individuals), event[2])[2]
    rates.infection.external[individual] = 0.
    rates.infection.internal[:, individual] = 0.
    @simd for i in find(states.susceptible)
      rates.infection.internal[individual, i] = riskfuncs.susceptibility(riskparams.susceptibility, population, i) *
                                                riskfuncs.transmissibility(riskparams.transmissibility, population, individual) *
                                                riskfuncs.infectivity(riskparams.infectivity, population, i, individual)
    end
  end
  return rates
end
