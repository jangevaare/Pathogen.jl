"""
loglikelihood(network::Network,
              networkrates::NetworkRates)
"""
function loglikelihood(network::Network,
                       networkrates::NetworkRates)
  # Initialize
  ll = 0.
  exposed = find([any(network.internal[:, i]) | network.external[i] for i = 1:length(network.external)])
  for i in exposed
    ll == -Inf && break
    # Numerator
    ll = log(networkrates.external[network.external[i]] + sum(networkrates.internal[:, network.internal[:, i]]))
    # Denominator
    ll -= log(networkrates.external[i] + sum(networkrates.internal[:, i]))
  end
  return ll
end


"""
loglikelihood(riskparams::SEIR_RiskParameters,
              events::SEIR_Events,
              riskfuncs::SEIR_RiskFunctions,
              population::DataFrame)

Calculates the log likelihood of a continuous time individual level model of
infectious disease transmission
"""
function loglikelihood(riskparams::SEIR_RiskParameters,
                       events::SEIR_Events,
                       riskfuncs::SEIR_RiskFunctions,
                       population::DataFrame)
  # Initialize
  ll = 0.
  eventtimes = [events.exposed events.infected events.removed]
  states = SEIR_States(population)
  rates = initialize_rates(states, population, riskfuncs, riskparams)
  network_rates = NetworkRates(events.individuals)

  # Find event order
  eventorder = sortperm(eventtimes[:])

  for i = 1:length(eventorder)
    # Stop log likelihood calculation after the last event
    isnan(eventtimes[eventorder[i]]) && break

    # Stop log likelihood calculation anytime the loglikelihood goes to -Inf
    ll == -Inf && break

    # Convert linear index to an event tuple (individual, event type)
    individual, eventtype = ind2sub(size(eventtimes), eventorder[i])

    # Find the rate total
    ratetotal = sum([sum(rates.exposure.external);
                     sum(rates.exposure.internal);
                     sum(rates.infection);
                     sum(rates.removal)])

    if i > 1
      # Find the time difference between consecutive events
      ΔT = eventtimes[eventorder[i]] - eventtimes[eventorder[i-1]]

      # loglikelihood contribution of specific event
      ll += log(ratetotal) - ratetotal * ΔT
    end

    # For exposure events
    if eventtype == 1
      # Copy exposure rates from this moment in time
      network_rates.external[individual] = copy(rates.exposure.external[individual])
      network_rates.internal[:, individual] = copy(rates.exposure.internal[:, individual])
      exposuretotal = network_rates.external[individual] + sum(network_rates.internal[:, individual])

      # loglikelihood contribution of an exposure event
      ll += log(exposuretotal/ratetotal)
      update_states!(states, SEIR_Event(1, individual))
      update_rates!(rates, states, SEIR_Event(1, individual), population, riskfuncs, riskparams)
    # For non-exposure events
    else
      # loglikelihood contribution of a non-exposure event
      ll += log(rates[eventtype][individual]/ratetotal)
      update_states!(states, SEIR_Event(eventtype + 1, individual))
      update_rates!(rates, states, SEIR_Event(eventtype + 1, individual), population, riskfuncs, riskparams)
    end
  end
  return ll, network_rates
end


"""
loglikelihood(riskparams::SIR_RiskParameters,
              events::SIR_Events,
              riskfuncs::SIR_RiskFunctions,
              population::DataFrame)

Calculates the log likelihood of a continuous time individual level model of
infectious disease transmission
"""
function loglikelihood(riskparams::SIR_RiskParameters,
                       events::SIR_Events,
                       riskfuncs::SIR_RiskFunctions,
                       population::DataFrame)
  # Initialize
  ll = 0.
  eventtimes = [events.infected events.removed]
  states = SIR_States(population)
  rates = initialize_rates(states, population, riskfuncs, riskparams)
  network_rates = NetworkRates(events.individuals)

  # Find event order
  eventorder = sortperm(eventtimes[:])

  for i = 1:length(eventorder)
    # Stop log likelihood calculation after the last event
    isnan(eventtimes[eventorder[i]]) && break

    # Stop log likelihood calculation anytime the loglikelihood goes to -Inf
    ll == -Inf && break

    # Convert linear index to an event tuple (individual, event type)
    individual, eventtype = ind2sub(size(eventtimes), eventorder[i])

    # Find the rate total
    ratetotal = sum([sum(rates.infection.external);
                     sum(rates.infection.internal);
                     sum(rates.removal)])

    if i > 1
      # Find the time difference between consecutive events
      ΔT = eventtimes[eventorder[i]] - eventtimes[eventorder[i-1]]

      # loglikelihood contribution of specific event
      ll += log(ratetotal) - ratetotal * ΔT
    end

    # For infection events
    if eventtype == 1
      # Copy exposure rates from this moment in time
      network_rates.external[individual] = copy(rates.infection.external[individual])
      network_rates.internal[:, individual] = copy(rates.infection.internal[:, individual])
      exposuretotal = network_rates.external[individual] + sum(network_rates.internal[:, individual])

      # loglikelihood contribution of an exposure event
      ll += log(exposuretotal/ratetotal)
      update_states!(states, SIR_Event(1, individual))
      update_rates!(rates, states, SIR_Event(1, individual), population, riskfuncs, riskparams)
    # For non-exposure events
    else
      # loglikelihood contribution of a non-exposure event
      ll += log(rates[eventtype][individual]/ratetotal)
      update_states!(states, SIR_Event(eventtype+1, individual))
      update_rates!(rates, states, SIR_Event(eventtype+1, individual), population, riskfuncs, riskparams)
    end
  end
  return ll, network_rates
end


"""
loglikelihood(riskparams::SEI_RiskParameters,
              events::SEI_Events,
              riskfuncs::SEI_RiskFunctions,
              population::DataFrame)

Calculates the log likelihood of a continuous time individual level model of
infectious disease transmission
"""
function loglikelihood(riskparams::SEI_RiskParameters,
                       events::SEI_Events,
                       riskfuncs::SEI_RiskFunctions,
                       population::DataFrame)
  # Initialize
  ll = 0.
  eventtimes = [events.exposed events.infected]
  states = SEI_States(population)
  rates = initialize_rates(states, population, riskfuncs, riskparams)
  network_rates = NetworkRates(events.individuals)

  # Find event order
  eventorder = sortperm(eventtimes[:])

  for i = 1:length(eventorder)
    # Stop log likelihood calculation after the last event
    isnan(eventtimes[eventorder[i]]) && break

    # Stop log likelihood calculation anytime the loglikelihood goes to -Inf
    ll == -Inf && break

    # Convert linear index to an event tuple (individual, event type)
    individual, eventtype = ind2sub(size(eventtimes), eventorder[i])

    # Find the rate total
    ratetotal = sum([sum(rates.exposure.external);
                     sum(rates.exposure.internal);
                     sum(rates.infection)])

    if i > 1
      # Find the time difference between consecutive events
      ΔT = eventtimes[eventorder[i]] - eventtimes[eventorder[i-1]]

      # loglikelihood contribution of specific event
      ll += log(ratetotal) - ratetotal * ΔT
    end

    # For exposure events
    if eventtype == 1
      # Copy exposure rates from this moment in time
      network_rates.external[individual] = copy(rates.exposure.external[individual])
      network_rates.internal[:, individual] = copy(rates.exposure.internal[:, individual])
      exposuretotal = network_rates.external[individual] + sum(network_rates.internal[:, individual])

      # loglikelihood contribution of an exposure event
      ll += log(exposuretotal/ratetotal)
      update_states!(states, SEI_Event(1, individual))
      update_rates!(rates, states, SEI_Event(1, individual), population, riskfuncs, riskparams)
    # For non-exposure events
    else
      # loglikelihood contribution of a non-exposure event
      ll += log(rates[eventtype][individual]/ratetotal)
      update_states!(states, SEI_Event(eventtype+1, individual))
      update_rates!(rates, states, SEI_Event(eventtype+1, individual), population, riskfuncs, riskparams)
    end
  end
  return ll, network_rates
end


"""
loglikelihood(riskparams::SI_RiskParameters,
              events::SI_Events,
              riskfuncs::SI_RiskFunctions,
              population::DataFrame)

Calculates the log likelihood of a continuous time individual level model of
infectious disease transmission
"""
function loglikelihood(riskparams::SI_RiskParameters,
                       events::SI_Events,
                       riskfuncs::SI_RiskFunctions,
                       population::DataFrame)
  # Initialize
  ll = 0.
  eventtimes = events.infected
  states = SI_States(population)
  rates = initialize_rates(states, population, riskfuncs, riskparams)
  network_rates = NetworkRates(events.individuals)

  # Find event order
  eventorder = sortperm(eventtimes[:])

  for i = 1:length(eventorder)
    # Stop log likelihood calculation after the last event
    isnan(eventtimes[eventorder[i]]) && break

    # Stop log likelihood calculation anytime the loglikelihood goes to -Inf
    ll == -Inf && break

    # Convert linear index to an event tuple (individual, event type)
    individual = eventorder[i]

    # Find the rate total
    ratetotal = sum([sum(rates.infection.external);
                     sum(rates.infection.internal)])

    if i > 1
      # Find the time difference between consecutive events
      ΔT = eventtimes[eventorder[i]] - eventtimes[eventorder[i-1]]

      # loglikelihood contribution of specific event
      ll += log(ratetotal) - ratetotal * ΔT
    end

    # For infection events
    # Copy exposure rates from this moment in time
    network_rates.external[individual] = copy(rates.infection.external[individual])
    network_rates.internal[:, individual] = copy(rates.infection.internal[:, individual])
    exposuretotal = network_rates.external[individual] + sum(network_rates.internal[:, individual])

    # loglikelihood contribution of an exposure event
    ll += log(exposuretotal/ratetotal)
    update_states!(states, SI_Event(1, individual))
    update_rates!(rates, states, SI_Event(1, individual), population, riskfuncs, riskparams)
  end
  return ll, network_rates
end
