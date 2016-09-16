function loglikelihood(iter::PathogenIteration
                       riskfuncs::RiskFunctions,
                       population::DataFrame)
  # Initialize
  ll = 0.
  eventtimes = [iter.events.exposed
                iter.events.infected
                iter.events.removed]
  riskparams = iter.risk_parameters
  rates = initialize_rates(population, riskfuncs, riskparms)
  networkrates = [fill(0., size(population, 1)), fill(0., (size(population, 1), size(population, 1)))]

  # Find event order
  eventorder = sortperm(events[:])

  # Sum the log likelihood of each event, taking histories into account
  for i = 2:length(eventorder)

    # Stop log likelihood calculation after the last event
    isnan(eventtimes[eventorder[i]]) && break

    # Stop log likelihood calculation anytime the loglikelihood goes to -Inf
    if isnan(ll)
      ll = -Inf
    end
    ll == -Inf && break

    # Convert linear index to an event tuple (individual, event type)
    individual, eventtype = ind2sub(size(eventtimes), eventorder[i])

    # Find the "master rate" through the sum of rate_array
    totals = [sum(rates[1]);
              sum(rates[2]);
              sum(rates[3]);
              sum(rates[4])]

    total = sum(totals)

    # Find the time difference between consecutive events
    deltaT = eventimes[eventorder[i]] - eventimes[eventorder[i-1]]

    # loglikelihood contribution of event time
    ll += loglikelihood(Exponential(1/total), deltaT)

    # loglikelihood contribution of specific event
    if eventtype == 1
      networkrates[1][individual] = rates[1][individual]
      networkrates[2][:, individual] = rates[2][:, individual]
      exposuretotal = rates[1][individual] + sum(rates[2][:, individual])
      ll += log(exposuretotal/total)
      # Use Gibbs step to generate network, thus next line is commented
      # ll += log((networkrates[1][network[id[1]]] + networkrates[2][network[:, id[1]]])/exposuretotal)
      if network.interal[individual]
        update_rates(rates, (1, individual))
      else
        source = findfirst(network.external[:,individual])
        update_rates(rates, (2, sub2ind(size(network.external), source, individual)))
      end
      update_rates!(rates, )
    else
      ll += log(rates[eventtype+1][individual]/total)
      update_rates!(rates, (eventtype+1, individual))
    end
  end
  return ll, networkrates
end
