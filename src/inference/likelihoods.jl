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

    # Find the rate total
    total = sum([sum(rates[1]);
                 sum(rates[2]);
                 sum(rates[3]);
                 sum(rates[4])])

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
      if network.internal[individual]
        update_rates(rates, (1, individual))
      else
        source = findfirst(network.external[:,individual])
        update_rates(rates, (2, sub2ind(size(network.external), source, individual)))
      end
    else
      ll += log(rates[eventtype+1][individual]/total)
      update_rates!(rates, (eventtype+1, individual))
    end
  end
  return ll, networkrates
end
