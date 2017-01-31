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


function show(io::IO, object::PathogenTrace)
  print(io, "PathogenTrace object (MCMC iterations: $(length(object)))")
end


function show(io::IO, object::PathogenIteration)
  print(io, "MCMC iteration object")
end


"""
Metropolis-Hastings rejection using log posteriors
"""
function MHreject(lp1::Float64, lp2::Float64)
  return rand() > exp(lp1 - lp2)
end


"""
Metropolis-Hastings acceptance using log posteriors
"""
function MHaccept(lp1::Float64, lp2::Float64)
  return rand() <= exp(lp1 - lp2)
end

"""
Initialize MCMC
"""
function initialize_mcmc(event_obs::EventObservations,
                         seq_obs::Vector{Sequence},
                         events_proposal::Events
                         riskparameter_priors::RiskParameterPriors,
                         riskfuncs::RiskFunctions,
                         substitutionmodel_priors::SubstitutionModelPrior,
                         population::DataFrame)
  riskparameter_proposal = rand(riskparameter_priors)
  substitutionmodel_proposal = rand(substitutionmodel_priors)

  lprior = logprior(riskparameter_priors,
                    riskparameter_proposal)
  lprior += logprior(substitutionmodel_priors,
                     substitutionmodel_proposal)

  llikelihood, network_rates = loglikelihood(riskparameter_proposal,
                                             events_proposal,
                                             riskfuncs,
                                             population)
  network_proposal = rand(network_rates)
  tree_proposal = Tree{Sequence, Void}()
  generatetree!(tree_proposal,
                events_proposal,
                event_obs,
                network_proposal)
  leaves = findleaves(tree_proposal)
  for j = 1:length(leaves)
    setdata!(tree_proposal.nodes[leaves[j]], seq_obs[j])
  end
  llikelihood += loglikelihood(tree_proposal,
                               substitutionmodel_proposal)
  lposterior = lprior + llikelihood
  pathogen_trace = PathogenTrace([riskparameter_proposal],
                                 [events_proposal],
                                 [network_proposal],
                                 [lposterior])
  phylo_trace = PhyloTrace([substitutionmodel_proposal],
                           [tree_proposal],
                           [lposterior])
  return pathogen_trace, phylo_trace
end


"""
Phylodynamic ILM MCMC
"""
function mcmc!(pathogen_trace::PathogenTrace,
               phylo_trace::PhyloTrace,
               n::Int64,
               ILM_kernel_variance::Array{Float64, 2},
               phylogenetic_kernel_variance::Array{Float64, 2},
               event_variance::Float64,
               event_obs::EventObservations,
               seq_obs::Vector{Sequence},
               riskparameter_priors::RiskParameterPriors,
               riskfuncs::RiskFunctions,
               substitutionmodel_priors::SubstitutionModelPrior,
               population::DataFrame,
               iterprob::Vector{Float64})
  progressbar = Progress(n, 5, "Performing $n MCMC iterations...", 25)
  individuals = size(population, 1)
  acceptances = fill(false, (1, length(iterprob)))
  for i = 1:n
    next!(progressbar, showvalues = [("Iterations", size(acceptances, 1)); ("Acceptances", sum(acceptances, 1))])

    iterationtype = findfirst(rand(Multinomial(1, iterprob)))

    if iterationtype == 1
      riskparameter_proposal = propose(pathogen_trace.riskparameters[end],
                                       ILM_kernel_variance)
    else
      riskparameter_proposal = RiskParameters(pathogen_trace.riskparameters[end].sparks,
                                              pathogen_trace.riskparameters[end].susceptibility,
                                              pathogen_trace.riskparameters[end].transmissibility,
                                              pathogen_trace.riskparameters[end].infectivity,
                                              pathogen_trace.riskparameters[end].latency,
                                              pathogen_trace.riskparameters[end].removal)
    end

    if iterationtype == 2
      substitutionmodel_proposal = propose(phylo_trace.substitutionmodel[end],
                                           substitutionmodel_priors,
                                           phylogenetic_kernel_variance)
    else
      substitutionmodel_proposal = phylo_trace.substitutionmodel[end]
    end

    if iterationtype == 3
      events_proposal = propose(pathogen_trace.events[end],
                                pathogen_trace.network[end],
                                event_variance)
    else
      events_proposal = Events(pathogen_trace.events[end].exposed,
                               pathogen_trace.events[end].infected,
                               pathogen_trace.events[end].removed)
    end
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
      if iterationtype == 4
        network_proposal = rand(network_rates)
      else
        network_proposal = Network(pathogen_trace.network[end].external,
                                   pathogen_trace.network[end].internal)
      end
      tree_proposal = Tree{Sequence, Void}()
      generatetree!(tree_proposal,
                    events_proposal,
                    event_obs,
                    network_proposal)
      leaves = findleaves(tree_proposal)
      for j = 1:length(leaves)
        setdata!(tree_proposal.nodes[leaves[j]], seq_obs[j])
      end
      llikelihood += loglikelihood(tree_proposal,
                                   substitutionmodel_proposal)
    else
      llikelihood = -Inf
    end

    lposterior = lprior + llikelihood

    if MHaccept(lposterior, pathogen_trace.logposterior[end])
      push!(pathogen_trace, PathogenIteration(riskparameter_proposal,
                                              events_proposal,
                                              network_proposal,
                                              lposterior))
      push!(phylo_trace, PhyloIteration(substitutionmodel_proposal,
                                        tree_proposal,
                                        lposterior))
      acceptance = fill(false, (1, length(iterprob)))
      acceptance[iterationtype] = true
      acceptances = vcat(acceptances, acceptance)
    else
      push!(pathogen_trace, PathogenIteration(RiskParameters(pathogen_trace.riskparameters[end].sparks,
                                                             pathogen_trace.riskparameters[end].susceptibility,
                                                             pathogen_trace.riskparameters[end].transmissibility,
                                                             pathogen_trace.riskparameters[end].infectivity,
                                                             pathogen_trace.riskparameters[end].latency,
                                                             pathogen_trace.riskparameters[end].removal),
                                              Events(pathogen_trace.events[end].exposed,
                                                     pathogen_trace.events[end].infected,
                                                     pathogen_trace.events[end].removed),
                                              Network(pathogen_trace.network[end].external,
                                                      pathogen_trace.network[end].internal),
                                              pathogen_trace.logposterior[end]))
      push!(phylo_trace, PhyloIteration(phylo_trace.substitutionmodel[end],
                                        phylo_trace.tree[end],
                                        phylo_trace.logposterior[end]))
      acceptance = fill(false, (1, length(iterprob)))
      acceptances = vcat(acceptances, acceptance)
    end
  end
  return pathogen_trace, phylo_trace
end
