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
  riskparameter_proposal = rand(riskparameter_priors)
  substitutionmodel_proposal = rand(substitutionmodel_priors)
  events_proposal = rand(event_priors)
  lprior = logprior(riskparameter_priors,
                    riskparameter_proposal)
  lprior += logprior(substitutionmodel_priors,
                     substitutionmodel_proposal)
  lprior += logprior(event_priors,
                     events_proposal)
  llikelihood, network_rates1 = loglikelihood(riskparameter_proposal,
                                              events_proposal,
                                              riskfuncs,
                                              population)
  network_proposal1 = rand(network_rates)
  tree_proposal = generatetree(events_proposal,
                                event_obs,
                                network_proposal)
  llikelihood += loglikelihood(seq_obs,
                               tree_proposal,
                               substitutionmodel_proposal)
  lposterior = lprior + llikelihood
  pathogen_trace = PathogenTrace([riskparameter_proposal],
                                 [events_proposal],
                                 [network_proposal],
                                 [lposterior])
  phylo_trace = PhyloTrace([substitutionmodel_proposal],
                           [tree_proposal],
                           [lposterior])

  for i = 2:n
    next!(progressbar)
    riskparameter_proposal = propose(pathogen_trace.riskparameters[end],
                                     riskparameter_priors,
                                     ILM_kernel_variance)
    substitutionmodel_proposal = propose(phylo_trace.substitutionmodel[end],
                                         substitutionmodel_priors,
                                         phylogenetic_kernel_variance)
    updated_inds = sample(1:individuals, rand(Poisson(1.)))
    events_proposal = propose(updated_inds,
                              pathogen_trace.events[end],
                              event_priors,
                              pathogen_trace.network[end])
    # events_proposal = rand(event_priors)
    lprior = logprior(riskparameter_priors,
                      riskparameter_proposal)
    lprior += logprior(substitutionmodel_priors,
                       substitutionmodel_proposal)
    lprior += logprior(event_priors,
                       events_proposal)
    if lprior > -Inf
      llikelihood, network_rates = loglikelihood(riskparameter_proposal,
                                                 events_proposal,
                                                 riskfuncs,
                                                 population)
      # updated_inds = sample(1:individuals, rand(Poisson(3.)))
      # network_proposal1 = propose(updated_inds,
      #                             network_proposal2,
      #                             network_rates1)
      network_proposal = rand(network_rates)
      tree_proposal = generatetree(events_proposal,
                                   event_obs,
                                   network_proposal)
      llikelihood += loglikelihood(seq_obs,
                                   tree_proposal,
                                   substitutionmodel_proposal)
    end
    lposterior = lprior + llikelihood
    if MHreject(lposterior, pathogen_trace.logposterior[end])
      push!(pathogen_trace, PathogenIteration(pathogen_trace.riskparameters[end],
                                              pathogen_trace.events[end],
                                              pathogen_trace.network[end],
                                              pathogen_trace.logposterior[end]))
      push!(phylo_trace, PhyloIteration(phylo_trace.substitutionmodel[end],
                                        phylo_trace.tree[end],
                                        phylo_trace.logposterior[end]))
    else
      push!(pathogen_trace, PathogenIteration(riskparameter_proposal,
                                              events_proposal,
                                              network_proposal,
                                              lposterior))
      push!(phylo_trace, PhyloIteration(substitutionmodel_proposal,
                                        tree_proposal,
                                        lposterior))
    end
  end
  return pathogen_trace, phylo_trace
end
