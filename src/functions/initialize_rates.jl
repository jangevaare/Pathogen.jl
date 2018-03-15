"""
initialize_rates(states::SEIR_States,
                 population::DataFrame,
                 riskfuncs::SEIR_RiskFunctions,
                 riskparams::SEIR_RiskParameters)

Create and initialize a rate array
"""
function initialize_rates(states::SEIR_States,
                          population::DataFrame,
                          riskfuncs::SEIR_RiskFunctions,
                          riskparams::SEIR_RiskParameters)
  individuals = size(population, 1)
  rates = SEIR_Rates(individuals)

  # Exposure
  for i in find(states.susceptible)
    # External exposure
    rates.exposure.external[i] = riskfuncs.susceptibility(riskparams.susceptibility, population, i) * riskfuncs.sparks(riskparams.sparks, population, i)

    # Internal exposure
    @simd for k in find(states.infected)
      rates.exposure.internal[k, i] = riskfuncs.susceptibility(riskparams.susceptibility, population, i) *
                                      riskfuncs.transmissibility(riskparams.transmissibility, population, k) *
                                      riskfuncs.infectivity(riskparams.infectivity, population, i, k)
    end
  end

  # Infection onset
  @simd for j in find(states.exposed)
    rates.infection[j] = riskfuncs.latency(riskparams.latency, population, j)
  end

  # Removal
  @simd for k in find(states.infected)
    rates.removal = riskfuncs.removal(riskparams.removal, population, k)
  end

  return rates
end


"""
initialize_rates(states::SIR_States,
                 population::DataFrame,
                 riskfuncs::SIR_RiskFunctions,
                 riskparams::SIR_RiskParameters)

Create and initialize a rate array
"""
function initialize_rates(states::SIR_States,
                          population::DataFrame,
                          riskfuncs::SIR_RiskFunctions,
                          riskparams::SIR_RiskParameters)
  individuals = size(population, 1)
  rates = SIR_Rates(individuals)

  # Infection
  for i in find(states.susceptible)
    # External infection
    rates.infection.external[i] = riskfuncs.susceptibility(riskparams.susceptibility, population, i) * riskfuncs.sparks(riskparams.sparks, population, i)

    # Internal infection
    @simd for k in find(states.infected)
      rates.infection.internal[k, i] = riskfuncs.susceptibility(riskparams.susceptibility, population, i) *
                                       riskfuncs.transmissibility(riskparams.transmissibility, population, k) *
                                       riskfuncs.infectivity(riskparams.infectivity, population, i, k)
    end
  end

  # Removal
  @simd for k in find(states.infected)
    rates.removal = riskfuncs.removal(riskparams.removal, population, k)
  end

  return rates
end


"""
initialize_rates(states::SEI_States,
                 population::DataFrame,
                 riskfuncs::SEI_RiskFunctions,
                 riskparams::SEI_RiskParameters)

Create and initialize a rate array
"""
function initialize_rates(states::SEI_States,
                          population::DataFrame,
                          riskfuncs::SEI_RiskFunctions,
                          riskparams::SEI_RiskParameters)
  individuals = size(population, 1)
  rates = SEI_Rates(individuals)

  # Exposure
  for i in find(states.susceptible)
    # External exposure
    rates.exposure.external[i] = riskfuncs.susceptibility(riskparams.susceptibility, population, i) * riskfuncs.sparks(riskparams.sparks, population, i)

    # Internal exposure
    @simd for k in find(states.infected)
      rates.exposure.internal[k, i] = riskfuncs.susceptibility(riskparams.susceptibility, population, i) *
                                      riskfuncs.transmissibility(riskparams.transmissibility, population, k) *
                                      riskfuncs.infectivity(riskparams.infectivity, population, i, k)
    end
  end

  # Infection onset
  @simd for j in find(states.exposed)
    rates.infection[j] = riskfuncs.latency(riskparams.latency, population, j)
  end

  return rates
end


"""
initialize_rates(states::SI_States,
                 population::DataFrame,
                 riskfuncs::SI_RiskFunctions,
                 riskparams::SI_RiskParameters)

Create and initialize a rate array
"""
function initialize_rates(states::SI_States,
                          population::DataFrame,
                          riskfuncs::SI_RiskFunctions,
                          riskparams::SI_RiskParameters)
  individuals = size(population, 1)
  rates = SI_Rates(individuals)

  # Infection
  for i in find(states.susceptible)
    # External infection
    rates.infection.external[i] = riskfuncs.susceptibility(riskparams.susceptibility, population, i) * riskfuncs.sparks(riskparams.sparks, population, i)

    # Internal infection
    @simd for k in find(states.infected)
      rates.infection.internal[k, i] = riskfuncs.susceptibility(riskparams.susceptibility, population, i) *
                                       riskfuncs.transmissibility(riskparams.transmissibility, population, k) *
                                       riskfuncs.infectivity(riskparams.infectivity, population, i, k)
    end
  end
  return rates
end
