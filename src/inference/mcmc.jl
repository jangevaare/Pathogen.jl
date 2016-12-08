"""
MCMC iteration
"""
type PathogenIteration
  riskparameters::RiskParameters
  events::Events
  network::Network
  logposterior::Float64
  acceptance::Bool
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
  acceptance::Vector{Bool}
end


function push!(trace::PathogenTrace, iteration::PathogenIteration)
  push!(trace.riskparameters, iteration.riskparameters)
  push!(trace.events, iteration.events)
  push!(trace.network, iteration.network)
  push!(trace.logposterior, iteration.logposterior)
  push!(trace.acceptance, iteration.acceptance)
  return trace
end


function append!(trace1::PathogenTrace, trace2::PathogenTrace)
  append!(trace1.riskparameters, trace2.riskparameters)
  append!(trace1.events, trace2.events)
  append!(trace1.network, trace2.network)
  append!(trace1.logposterior, trace2.logposterior)
  append!(trace1.acceptance, trace2.acceptance)
  return trace1
end


function length(x::PathogenTrace)
  return length(x.riskparameters)
end


function show(io::IO, object::PathogenTrace)
  print(io, "PathogenTrace object (MCMC iterations: $(length(object)), acceptance rate: $(trunc(sum(object.acceptance)*100/length(object), 4))%)")
end


function show(io::IO, object::PathogenIteration)
  print(io, "MCMC iteration object")
end


"""
Metropolis-Hastings rejection using log posteriors
"""
function MHreject(lp1::Float64, lp2::Float64)
  if rand() <= exp(lp1 - lp2)
    return false
  else
    return true
  end
end


"""
Metropolis-Hastings acceptance using log posteriors
"""
function MHaccept(lp1::Float64, lp2::Float64)
  if rand() < exp(lp1 - lp2)
    return true
  else
    return false
  end
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
  llikelihood, network_rates = loglikelihood(riskparameter_proposal,
                                             events_proposal,
                                             riskfuncs,
                                             population)
  network_proposal = rand(network_rates)
  tree_proposal = generatetree(events_proposal,
                               event_obs,
                               network_proposal)
  leaves = findleaves(tree_proposal)
  for j = 1:length(leaves)
    setdata!(tree_proposal.nodes[leave[j]], seq_obs[j])
  end
  llikelihood += loglikelihood(tree_proposal,
                              substitutionmodel_proposal)
  lposterior = lprior + llikelihood
  pathogen_trace = PathogenTrace([riskparameter_proposal],
                                 [events_proposal],
                                 [network_proposal],
                                 [lposterior],
                                 [true])
  phylo_trace = PhyloTrace([substitutionmodel_proposal],
                           [tree_proposal],
                           [lposterior],
                           [true])

  for i = 2:n
    next!(progressbar)
    # iterationtype = findfirst(rand(Multinomial(1, [0.25; 0.25; 0.25; 0.25])))
    iterationtype = findfirst(rand(Multinomial(1, [1.0; 0.0; 0.0; 0.0])))

    if iterationtype == 1
      riskparameter_proposal = propose(pathogen_trace.riskparameters[end],
                                       riskparameter_priors,
                                       ILM_kernel_variance)
    else
      riskparameter_proposal = pathogen_trace.riskparameters[end]
    end

    if iterationtype == 2
      substitutionmodel_proposal = propose(phylo_trace.substitutionmodel[end],
                                           substitutionmodel_priors,
                                           phylogenetic_kernel_variance)
    else
      substitutionmodel_proposal = phylo_trace.substitutionmodel[end]
    end

    if iterationtype == 3
      updated_inds = sample(1:individuals, 1)
      events_proposal = propose(updated_inds,
                                pathogen_trace.events[end],
                                event_priors,
                                pathogen_trace.network[end])
      # events_proposal = rand(event_priors)
    else
      events_proposal = pathogen_trace.events[end]
    end

    lprior = logprior(riskparameter_priors,
                      riskparameter_proposal)
    lprior += logprior(substitutionmodel_priors,
                       substitutionmodel_proposal)
    lprior += logprior(event_priors,
                       events_proposal)

    llikelihood, network_rates = loglikelihood(riskparameter_proposal,
                                               events_proposal,
                                               riskfuncs,
                                               population)
    if iterationtype == 4
      updated_inds = sample(1:individuals, 1)
      network_proposal = propose(updated_inds,
                                 pathogen_trace.network[end],
                                 network_rates)
      # network_proposal = rand(network_rates)
    else
      network_proposal = pathogen_trace.network[end]
    end

    tree_proposal = generatetree(events_proposal,
                                 event_obs,
                                 network_proposal)
    leaves = findleaves(tree_proposal)
    for j = 1:length(leaves)
      setdata!(tree_proposal.nodes[leave[j]], seq_obs[j])
    end
    llikelihood += loglikelihood(tree_proposal,
                                 substitutionmodel_proposal)

    lposterior = lprior + llikelihood

    if MHaccept(lposterior, pathogen_trace.logposterior[end])
      push!(pathogen_trace, PathogenIteration(riskparameter_proposal,
                                              events_proposal,
                                              network_proposal,
                                              lposterior,
                                              true))
      push!(phylo_trace, PhyloIteration(substitutionmodel_proposal,
                                        tree_proposal,
                                        lposterior,
                                        true))
    else
      push!(pathogen_trace, PathogenIteration(pathogen_trace.riskparameters[end],
                                              pathogen_trace.events[end],
                                              pathogen_trace.network[end],
                                              pathogen_trace.logposterior[end],
                                              false))
      push!(phylo_trace, PhyloIteration(phylo_trace.substitutionmodel[end],
                                        phylo_trace.tree[end],
                                        phylo_trace.logposterior[end],
                                        false))
    end
  end
  return pathogen_trace, phylo_trace
end
