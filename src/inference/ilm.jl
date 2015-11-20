"""
Calculate the log likelihood and return an exposure network array under specified parameters values and observations

α, β: powerlaw exposure kernel parameters
η: external pressure rate
ρ: infectivity rate (1/mean latent period)
γ: removal rate (1/mean infectious period)
"""
function SEIR_loglikelihood(α::Float64,
                            β::Float64,
                            η::Float64,
                            ρ::Float64,
                            γ::Float64,
                            aug::SEIR_augmented,
                            obs::SEIR_observed,
                            debug=false::Bool,
                            dist=Euclidean())

  # Initiate an exposure network
  network_rates = fill(0., (1 + length(obs.covariates), length(obs.covariates)))

  # Start with loglikelihood of 0
  ll = 0.

  # Create event timing array
  event_times = [aug.exposed aug.infectious aug.removed]

  # Find event order
  event_order = sortperm(event_times[:])

  # Create empty rate array
  rate_array = fill(0., (1 + length(obs.covariates) + 2, length(obs.covariates)))

  # First row is external pressure rate
  rate_array[1, :] = η

  # Sum the log likelihood of each event, taking histories into account
  for i = 1:length(event_order)

    # Stop log likelihood calculation after the last event
    isnan(event_times[event_order[i]]) && break

    # Stop log likelihood calculation anytime the loglikelihood goes to -Inf
    if isnan(ll)
      ll = -Inf
    end
    ll == -Inf && break

    # Convert linear index to an event tuple (individual, event type)
    id = ind2sub(size(event_times), event_order[i])

    # Find the "master rate" through the sum of rate_array
    master_rate = sum(rate_array)

    # Don't consider likelilihood contribution of first event time
    if i > 1
      # loglikelihood of event time with master rate
      ll += logpdf(Exponential(1.0/master_rate), event_times[event_order[i]] - event_times[event_order[i-1]])
    end

    # loglikelihood of the one event that did occur
    ll += log(sum(rate_array[:, id[1]])/master_rate)

    # if debug
    #   println("Probability of event $i: $(sum(rate_array[:, id[1]])/master_rate)")
    # end

    # Exposure event
    if id[2] == 1

      # Record exposure rates at time of exposure
      network_rates[:, id[1]] = rate_array[1:(size(rate_array, 2) + 1), id[1]]

      # Update exposure rates
      rate_array[1:(size(rate_array, 2) + 1), id[1]] = 0.

      # Update infectivity rate
      rate_array[1 + size(rate_array, 2) + 1, id[1]] = ρ

    # Infectiousness event
    elseif id[2] == 2

      # Update infectivity rate
      rate_array[1 + size(rate_array, 2) + 1, id[1]] = 0.

      # Update removal rate
      rate_array[1 + size(rate_array, 2) + 2, id[1]] = γ

      # Update exposure rates for rest of susceptible population
      for j = 1:size(rate_array, 2)
        if j != id[1] && rate_array[1, j] > 0.
          rate_array[id[1] + 1, j] = α*evaluate(dist,
                                                obs.covariates[id[1]],
                                                obs.covariates[j])^-β
        end
      end

    # Removal event
    elseif id[2] == 3

      # Update removal rate
      rate_array[1 + size(rate_array, 2) + 2, id[1]] = 0.

      # Update exposure rates
      rate_array[id[1] + 1, :] = 0.
    end

    # Provide loop position when loglikelihood goes to -Inf when debugging
    if debug && ll == -Inf
      if id[2] == 1
        println("Event $i (exposure of individual $(id[1])) caused log likelihood to go to -Inf")
      elseif id[2] == 2
        println("Event $i (infection of individual $(id[1])) caused log likelihood to go to -Inf")
      elseif id[2] == 3
        println("Event $i (removal of individual $(id[1])) caused log likelihood to go to -Inf")
      end
    end
  end
  debug && println("SEIR loglikelihood: $ll")
  return ll, network_rates
end
