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
  lpgposterior::Vector{Float64}
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
              event_obs::EventObservations,
              seq_obs::Array{Bool, 3},
              event_priors::EventPriors,
              riskparameter_priors::RiskParameterPriors,
              riskfuncs::RiskFunctions,
              substitutionmodel_func::Function,
              substitutionmodel_priors::SubstitutionModelPriors,
              population::DataFrame,
              tune=1000::Int64)
  progressbar = Progress(n, 5, "Performing $n MCMC iterations...", 25)
  individuals = size(population, 1)
  transition_kernel_var1 = transition_kernel_variance(riskparameter_priors)
  transition_kernel_var2 = transition_kernel_variance(substitutionmodel_priors)
  riskparameter_proposal1 = rand(riskparameter_priors)
  substitutionmodel_proposal1 = substitutionmodel(rand(substitutionmodel_priors))
  lprior1 = logprior(riskparameter_proposal1,
                     riskparameter_priors)
  lprior1 += logprior(substitutionmodel_proposal1,
                      substitutionmodel_priors)
  events_proposal1 = rand(event_priors)
  llikelihood1, network_rates1 = loglikelihood(riskparameter_proposal1,
                                               events_proposal1,
                                               riskfuncs,
                                               population)
  network_proposal1 = rand(network_rates1)
  tree_proposal1, observednodes1 = generatetree(events_proposal1,
                                                event_obs,
                                                network_proposal1)
  llikelihood1 += loglikelihood(seq_obs,
                                tree_proposal1,
                                substitutionmodel_proposal1)
  lposterior1 = lprior1 + llikelihood1
  pathogen_trace = PathogenTrace([riskparameter_proposal1],
                                 [events_proposal1],
                                 [network_proposal1],
                                 [lposterior1])
  phylo_trace = PhyloTrace(substitutionmodel_proposal1,
                           tree_proposal1,
                           lposterior1)

  for i = 1:n
    next!(progressbar)
    if mod(tune, i) == 0
      transition_kernel_var1 = transition_kernel_variance(pathogen_trace.riskparameters)
      transition_kernel_var2 = transition_kernel_variance(phylo_trace.Θ)
    end

    riskparameter_proposal2 = riskparameter_proposal1
    substitutionmodel_proposal2 = substitutionmodel_proposal1
    lprior2 = lprior1
    events_proposal2 = events_proposal1
    network_proposal2 = network_proposal1
    tree_proposal2 = tree_proposal1
    lposterior2 = lposterior1

    if mod(3, i) == 0
      lprior1 = -Inf
      while lprior1 == -Inf
        riskparameter_proposal1 = propose(riskparameter_proposal2,
                                          transition_kernel_var1)
        substitutionmodel_proposal1 = propose(substitutionmodel_proposal2,
                                              transition_kernel_var2)
        lprior1 = logprior(riskparameter_proposal1,
                           riskparameter_priors)
        lprior1 += logprior(substitutionmodel_proposal1,
                            substitutionmodel_priors)
      end
    end

    if mod(3, i) == 1
      updated_inds = sample(1:individuals, rand(Poisson(1.)))
      events_proposal1 = propose(updated_inds,
                                 events_proposal2,
                                 event_priors,
                                 network_proposal1)
    end

    if mod(3, i) != 2
      llikelihood1, network_rates1 = loglikelihood(riskparameter_proposal1,
                                                   events_proposal1,
                                                   network_proposal1,
                                                   riskfuncs,
                                                   population)
    end
    if mod(3, i) == 2
      updated_inds = sample(1:individuals, rand(Poisson(1.)))
      network_proposal1 = propose(updated_inds,
                                  network_rates1,
                                  network_proposal2)
    end
    tree_proposal1, observednodes1 = generatetree(events_proposal1,
                                                  event_obs,
                                                  network_proposal1)
    llikelihood1 += loglikelihood(seq_obs,
                                  tree_proposal1,
                                  substitutionmodel_proposal1)
    lposterior1 = lprior1 + llikelihood1
    if MHreject(lposterior1, lposterior2)
      lprior1 = lprior2
      llikelihood1 = llikelihood2
      riskparameter_proposal1 = riskparameter_proposal2
      events_proposal1 = events_proposal2
      network_proposal1 = networkproposal2
      tree_proposal1 = tree_proposal2
      lposterior1 = lposterior2
    end
    push!(pathogen_trace, PathogenProposal(riskparameter_proposal1,
                                           events_proposal1,
                                           network_proposal1,
                                           lposterior1))
    push!(phylo_trace, PhyloTrace(substitutionmodel_proposal1,
                                  tree_proposal1,
                                  lposterior1))
  end
  return pathogen_trace, phylo_trace
end
