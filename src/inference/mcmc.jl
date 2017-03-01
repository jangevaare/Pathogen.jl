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


function deleteat!(trace::PathogenTrace, inds)
  deleteat!(trace.riskparameters, inds)
  deleteat!(trace.events, inds)
  deleteat!(trace.network, inds)
  deleteat!(trace.logposterior, inds)
  return trace
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
                         seq_obs::Dict{Int64, Sequence},
                         events_proposal::Events,
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
  tree_proposal = generate_tree(events_proposal,
                               event_obs,
                               network_proposal)
  llikelihood += loglikelihood(tree_proposal,
                               substitutionmodel_proposal,
                               seq_obs)
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
mcmc!(pathogen_trace::PathogenTrace,
      phylo_trace::PhyloTrace,
      n::Int64,
      thin::Int64,
      ILM_kernel_variance::Array{Float64, 2},
      phylogenetic_kernel_variance::Array{Float64, 2},
      event_variance::Float64,
      event_obs::EventObservations,
      seq_obs::Dict{Int64, Sequence},
      riskparameter_priors::RiskParameterPriors,
      riskfuncs::RiskFunctions,
      substitutionmodel_priors::SubstitutionModelPrior,
      population::DataFrame,
      iterprob::Vector{Float64})

Phylodynamic ILM MCMC
"""
function mcmc!(pathogen_trace::PathogenTrace,
               phylo_trace::PhyloTrace,
               n::Int64,
               thin::Int64,
               ILM_kernel_variance::Array{Float64, 2},
               phylogenetic_kernel_variance::Array{Float64, 2},
               event_variance::Float64,
               event_obs::EventObservations,
               seq_obs::Dict{Int64, Sequence},
               riskparameter_priors::RiskParameterPriors,
               riskfuncs::RiskFunctions,
               substitutionmodel_priors::SubstitutionModelPrior,
               population::DataFrame,
               iterprob::Vector{Float64})
  if thin < 0
    error("Thinning rate can not be less than 1")
  end
  progressbar = Progress(n, 5, "Performing $n MCMC iterations...", 25)
  individuals = event_obs.individuals
  validevents = find([!isnan(event_obs.infected) !isnan(event_obs.infected) !isnan(event_obs.removed)])
  validexposures = find(!isnan(event_obs.infected))

  riskparameters_previous = copy(pathogen_trace.riskparameters[end])
  substitutionmodel_previous = copy(phylo_trace.substitutionmodel[end])
  events_previous = copy(pathogen_trace.events[end])
  network_previous = copy(pathogen_trace.network[end])
  tree_previous = phylo_trace.tree[end]
  lposterior_previous = copy(pathogen_trace.logposterior[end])

  for i = 1:n
    next!(progressbar)
    iterationtype = findfirst(rand(Multinomial(1, iterprob)))

    if iterationtype == 1
      riskparameter_proposal = propose(pathogen_trace.riskparameters[end],
                                       ILM_kernel_variance)
    else
      riskparameter_proposal = copy(riskparameters_previous)
    end

    if iterationtype == 2
      substitutionmodel_proposal = propose(phylo_trace.substitutionmodel[end],
                                           phylogenetic_kernel_variance)
    else
      substitutionmodel_proposal = copy(substitutionmodel_previous)
    end

    if iterationtype == 3
      i, j = ind2sub((individuals, 3), sample(validevents))
      events_proposal = propose(i, j,
                                events_previous,
                                network_previous,
                                event_obs,
                                event_variance)
    else
      events_proposal = copy(events_previous)
    end

    lprior = logprior(riskparameter_priors,
                      riskparameter_proposal)
    lprior += logprior(substitutionmodel_priors,
                       substitutionmodel_proposal)
    if lprior > -Inf
      llikelihood, network_rates = loglikelihood(riskparameter_proposal,
                                                 events_proposal,
                                                 riskfuncs,
                                                 population)
      if iterationtype == 4
        i = sample(validexposures)
        network_proposal = propose(i,
                                   network_previous,
                                   network_rates)
      elseif iterationtype == 5
        network_proposal = rand(network_rates)
      else
        network_proposal = copy(network_previous)
      end

      tree_proposal = generate_tree(events_proposal,
                                    event_obs,
                                    network_proposal)
      llikelihood += loglikelihood(tree_proposal,
                                   substitutionmodel_proposal,
                                   seq_obs)
    else
      llikelihood = -Inf
    end

    lposterior = lprior + llikelihood

    if MHaccept(lposterior, lposterior_previous)
      riskparameters_previous = riskparameter_proposal
      substitutionmodel_previous = substitutionmodel_proposal
      events_previous = events_proposal
      network_previous = network_proposal
      tree_previous = tree_proposal
      lposterior_previous = lposterior
    end
    if mod(i, thin) == 0
      push!(pathogen_trace, PathogenIteration(copy(riskparameters_previous),
                                              copy(events_previous),
                                              copy(network_previous),
                                              copy(lposterior_previous)))
      push!(phylo_trace, PhyloIteration(copy(substitutionmodel_previous),
                                        tree_previous,
                                        copy(lposterior_previous)))
    end
  end
  return pathogen_trace, phylo_trace
end
