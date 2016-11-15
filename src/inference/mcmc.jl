"""
MCMC iteration
"""
type PathogenIteration
  riskparameters::RiskParameters
  events::Events
  network::Network
  logposterior::Float64
end


PathogenProposal = PathogenIteration


"""
MCMC trace
"""
type PathogenTrace
  riskparameters::Vector{RiskParameters}
  events::Vector{Events}
  network::Vector{Network}
  logposterior::Vector{Float64}
end


function push!(trace::PathogenTrace, iteration::PathogenIteration)
  push!(trace.riskparameters, iteration.riskparameters)
  push!(trace.events, iteration.events)
  push!(trace.network, iteration.network)
  push!(trace.logposterior, iteration.logposterior)
  return trace
end


function append!(trace1::PathogenTrace, trace2::PathogenTrace)
  append!(trace1.riskparameters, trace2.riskparameters)
  append!(trace1.events, trace2.events)
  append!(trace1.network, trace2.network)
  append!(trace1.logposterior, trace2.logposterior)
  return trace1
end


function length(x::PathogenTrace)
  return length(x.riskparameters)
end


"""
Metropolis-Hastings rejection using log posteriors
"""
function MHreject(lp1::Float64, lp2::Float64)
  reject = true
  if lp1 >= lp2
    reject = false
  elseif exp(lp1 - lp2) >= rand()
    reject = false
  end
  return reject
end


"""
Metropolis-Hastings acceptance using log posteriors
"""
function MHaccept(lp1::Float64, lp2::Float64)
  return !MHreject(lp1, lp2)
end


"""
Phylodynamic ILM MCMC
"""
function mcmc(n::Int64,
              ILM_kernel_variance::Vector{Float64},
              phylogenetic_kernel_variance::Vector{Float64},
              event_obs::EventObservations,
              seq_obs::Vector{Sequence},
              event_priors::EventPriors,
              riskparameter_priors::RiskParameterPriors,
              riskfuncs::RiskFunctions,
              substitutionmodel_priors::SubstitutionModelPrior,
              population::DataFrame)
  progressbar = Progress(n, 5, "Performing $n MCMC iterations...", 25)
  individuals = size(population, 1)
  riskparameter_proposal1 = rand(riskparameter_priors)
  substitutionmodel_proposal1 = rand(substitutionmodel_priors)
  events_proposal1 = rand(event_priors)
  lprior = logprior(riskparameter_priors,
                    riskparameter_proposal1)
  lprior += logprior(substitutionmodel_priors,
                     substitutionmodel_proposal1)
  lprior += logprior(event_priors,
                     events_proposal1)
  llikelihood, network_rates1 = loglikelihood(riskparameter_proposal1,
                                              events_proposal1,
                                              riskfuncs,
                                              population)
  network_proposal1 = rand(network_rates1)
  tree_proposal1 = generatetree(events_proposal1,
                                event_obs,
                                network_proposal1)
  llikelihood += loglikelihood(seq_obs,
                               tree_proposal1,
                               substitutionmodel_proposal1)
  lposterior1 = lprior + llikelihood
  pathogen_trace = PathogenTrace([riskparameter_proposal1],
                                 [events_proposal1],
                                 [network_proposal1],
                                 [lposterior1])
  phylo_trace = PhyloTrace([substitutionmodel_proposal1],
                           [tree_proposal1],
                           [lposterior1])

  for i = 1:n
    next!(progressbar)
    riskparameter_proposal2 = riskparameter_proposal1
    substitutionmodel_proposal2 = substitutionmodel_proposal1
    events_proposal2 = events_proposal1
    network_proposal2 = network_proposal1
    tree_proposal2 = tree_proposal1
    lposterior2 = lposterior1

    if mod(i, 4) == 1
      riskparameter_proposal1 = propose(riskparameter_proposal2,
                                        riskparameter_priors,
                                        ILM_kernel_variance)
    end
    if mod(i, 4) == 2
      substitutionmodel_proposal1 = propose(substitutionmodel_proposal2,
                                            substitutionmodel_priors,
                                            phylogenetic_kernel_variance)
    end
    if mod(i, 4) == 3
      updated_inds = sample(1:individuals, rand(Poisson(1.)))
      events_proposal1 = propose(updated_inds,
                                 events_proposal2,
                                 event_priors,
                                 network_proposal1)
    end
    lprior = logprior(riskparameter_priors,
                      riskparameter_proposal1)
    lprior += logprior(substitutionmodel_priors,
                       substitutionmodel_proposal1)
    lprior += logprior(event_priors,
                       events_proposal1)
    if lprior > -Inf
      llikelihood, network_rates1 = loglikelihood(riskparameter_proposal1,
                                                  events_proposal1,
                                                  riskfuncs,
                                                  population)
      if mod(i, 4) == 0
        updated_inds = sample(1:individuals, rand(Poisson(1.)))
        network_proposal1 = propose(updated_inds,
                                    network_proposal2,
                                    network_rates1)
      end
      tree_proposal1 = generatetree(events_proposal1,
                                    event_obs,
                                    network_proposal1)
      llikelihood += loglikelihood(seq_obs,
                                   tree_proposal1,
                                   substitutionmodel_proposal1)
    end
    lposterior1 = lprior + llikelihood
    if MHreject(lposterior1, lposterior2)
      riskparameter_proposal1 = riskparameter_proposal2
      events_proposal1 = events_proposal2
      network_proposal1 = network_proposal2

      substitutionmodel_proposal1 = substitutionmodel_proposal2
      tree_proposal1 = tree_proposal2

      lposterior1 = lposterior2
    end
    push!(pathogen_trace, PathogenIteration(riskparameter_proposal1,
                                            events_proposal1,
                                            network_proposal1,
                                            lposterior1))
    push!(phylo_trace, PhyloIteration(substitutionmodel_proposal1,
                                      tree_proposal1,
                                      lposterior1))
  end
  return pathogen_trace, phylo_trace
end
