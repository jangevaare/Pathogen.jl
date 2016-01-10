"""
Propose a network
"""
function propose_network(changed_individuals::Vector{Int64},
                         network_rates::Array{Float64, 2},
                         previous_network::Array{Bool, 2},
                         debug=false::Bool)
  network = copy(previous_network)
  rate_totals = sum(network_rates,1)
  network[:, changed_individuals] = false
  @assert(size(network_rates) == size(previous_network),
          "A mismatch in the previous network and network rates dimensions was detected in the network proposal function")
  for i in changed_individuals
    network[findfirst(rand(Multinomial(1, network_rates[:,i]/rate_totals[i]))), i] = true
  end
  if debug
    println("Network proposal contains $(sum(network)) total exposure events, with up to $(length(changed_individuals)) change(s) from the previous network")
    println(spy(network))
  end
  return network
end


"""
Propose a network
"""
function propose_network(network_rates::Array{Float64, 2},
                         debug=false::Bool)
  network = fill(false, size(network_rates))
  rate_totals = sum(network_rates, 1)
  exposures = find(rate_totals .> 0)
  for i in exposures
    network[findfirst(rand(Multinomial(1, network_rates[:,i]/rate_totals[i]))), i] = true
  end
  if debug
    println("Network proposal ($(sum(network)) infections total):")
    println(spy(network))
  end
  return network
end


"""
Propose a network
"""
function propose_network(network_rates::Array{Float64, 2},
                         previous_network::Array{Bool, 2},
                         debug=false::Bool,
                         changes=rand(Poisson(1.))::Int64)
  if changes == 0
    return previous_network
  else
    rate_totals = sum(network_rates, 1)
    if changes >= sum(rate_totals .> 0)
      return propose_network(network_rates,
                             debug)
    else
      changed_individuals = sample(find(rate_totals .> 0), changes, replace=false)
      return propose_network(changed_individuals,
                             network_rates,
                             previous_network,
                             debug)
    end
  end
end


"""
For a given transmission network, find the time between the pathogen sequences between every individuals i and j
"""
function seq_distances(obs::SEIR_observed, aug::SEIR_augmented, network::Array{Bool, 2}, debug=false::Bool)
  pathways = pathwaysto(find(obs.sequenced), network, debug)
  seq_dist = fill(0., (size(network, 2), size(network, 2)))
  for i = 1:length(pathways)
    for j = 1:(i-1)
      k = 1
      while length(pathways[i]) > k && length(pathways[j]) > k && pathways[i][end - k] == pathways[j][end - k]
        k += 1
      end
      if k == length(pathways[i])
        seq_dist[pathways[i][1],pathways[j][1]] += obs.infectious[pathways[j][1]] - aug.exposed[pathways[j][end - k]]
        seq_dist[pathways[i][1],pathways[j][1]] += abs(aug.exposed[pathways[j][end - k]] - obs.infectious[pathways[i][1]])
      elseif k == length(pathways[j])
        seq_dist[pathways[i][1],pathways[j][1]] += obs.infectious[pathways[i][1]] - aug.exposed[pathways[i][end - k]]
        seq_dist[pathways[i][1],pathways[j][1]] += abs(aug.exposed[pathways[i][end - k]] - obs.infectious[pathways[j][1]])
      else
        seq_dist[pathways[i][1],pathways[j][1]] += obs.infectious[pathways[i][1]] - aug.exposed[pathways[i][end - k]]
        seq_dist[pathways[i][1],pathways[j][1]] += obs.infectious[pathways[j][1]] - aug.exposed[pathways[j][end - k]]
        seq_dist[pathways[i][1],pathways[j][1]] += abs(aug.exposed[pathways[j][end - k]] - aug.exposed[pathways[i][end - k]])
      end
    end
  end
  seq_dist += transpose(seq_dist)
  if debug
    @assert(all(seq_dist .>= 0), "Negative sequence distance detected...")
  end
  return seq_dist
end


"""
For a given transmission network, find the time between the pathogen sequences between every individuals i and j
"""
function seq_distances(obs::SIR_observed, aug::SIR_augmented, network::Array{Bool, 2}, debug=false::Bool)
  pathways = pathwaysto(find(obs.sequenced), network, debug)
  seq_dist = fill(0., (size(network, 2), size(network, 2)))
  for i = 1:length(pathways)
    for j = 1:(i-1)
      k = 1
      while length(pathways[i]) > k && length(pathways[j]) > k && pathways[i][end - k] == pathways[j][end - k]
        k += 1
      end
      if k == length(pathways[i])
        seq_dist[pathways[i][1],pathways[j][1]] += obs.infectious[pathways[j][1]] - aug.infectious[pathways[j][end - k]]
        seq_dist[pathways[i][1],pathways[j][1]] += abs(aug.infectious[pathways[j][end - k]] - obs.infectious[pathways[i][1]])
      elseif k == length(pathways[j])
        seq_dist[pathways[i][1],pathways[j][1]] += obs.infectious[pathways[i][1]] - aug.infectious[pathways[i][end - k]]
        seq_dist[pathways[i][1],pathways[j][1]] += abs(aug.infectious[pathways[i][end - k]] - obs.infectious[pathways[j][1]])
      else
        seq_dist[pathways[i][1],pathways[j][1]] += obs.infectious[pathways[i][1]] - aug.infectious[pathways[i][end - k]]
        seq_dist[pathways[i][1],pathways[j][1]] += obs.infectious[pathways[j][1]] - aug.infectious[pathways[j][end - k]]
        seq_dist[pathways[i][1],pathways[j][1]] += abs(aug.infectious[pathways[j][end - k]] - aug.infectious[pathways[i][end - k]])
      end
    end
  end
  seq_dist += transpose(seq_dist)
  if debug
    @assert(all(seq_dist .>= 0), "Negative sequence distance detected...")
  end
  return seq_dist
end


"""
Loglikelihood for a transmission network based on sequence data and event timing
"""
function phylogenetic_network_loglikelihood(obs::Observed,
                                            aug::Augmented,
                                            network::Array{Bool, 2},
                                            p_matrix::Function,
                                            debug=false::Bool)
  ll = 0.
  infected = find(obs.sequenced)
  seq_dist = seq_distances(obs, aug, network, debug)
  for i = 1:length(infected)
    for j = 1:(i-1)
      ll += sum(log(p_matrix(seq_dist[infected[i],infected[j]]))[sub2ind((4, 4), obs.seq[infected[i]], obs.seq[infected[j]])])
      (ll == -Inf || isnan(ll)) && break
    end
  end
  if isnan(ll)
    ll = -Inf
  end
  if debug
    println("Phylogenetic network log likelihood: $(round(ll,3))")
  end
  return ll
end


"""
Loglikelihood for a transmission network based on exposure rates
"""
function exposure_network_loglikelihood(network::Array{Bool, 2}, network_rates::Array{Float64}, debug=false::Bool)
  ll = 0.
  infected = find(sum(network, 1))
  for i in infected
     ll += log(network_rates[findfirst(network[:,i]),i]/sum(network_rates[:,i]))
     if debug
       ll == -Inf && println("Exposure of individual $i caused network -Inf log likelihood")
       isnan(ll) && println("Exposure of individual $i caused network NaN log likelihood")
     end
     (ll == -Inf || isnan(ll)) && break
  end
  if isnan(ll)
    ll = -Inf
  end
  if debug
    println("Exposure network log likelihood: $(round(ll, 3))")
  end
  return ll
end
