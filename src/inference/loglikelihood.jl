function loglikelihood(riskparams::RiskParameters,
                       events::Events,
                       riskfuncs::RiskFunctions,
                       population::DataFrame)
  # Initialize
  ll = 0.
  eventtimes = [events.exposed events.infected events.removed]
  rates = initialize_rates(population, riskfuncs, riskparams)
  networkrates = [fill(0., size(population, 1)), fill(0., (size(population, 1), size(population, 1)))]

  # Find event order
  eventorder = sortperm(eventtimes[:])

  for i = 1:length(eventorder)
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

    if i > 1
      # Find the time difference between consecutive events
      deltaT = eventtimes[eventorder[i]] - eventtimes[eventorder[i-1]]

      # loglikelihood contribution of event time
      ll += loglikelihood(Exponential(1/total), [deltaT])
    end

    # loglikelihood contribution of specific event
    if eventtype == 1
      networkrates[1][individual] = rates[1][individual]
      networkrates[2][:, individual] = rates[2][:, individual]
      exposuretotal = rates[1][individual] + sum(rates[2][:, individual])
      ll += log(exposuretotal/total)
      update_rates!(rates, (1, individual))
    else
      ll += log(rates[eventtype+1][individual]/total)
      update_rates!(rates, (eventtype+1, individual))
    end
  end
  return ll, networkrates
end
