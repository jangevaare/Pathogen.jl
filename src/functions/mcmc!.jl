"""
mcmc!(pathogen_trace::PathogenTrace,
      n::Int64,
      thin::Int64,
      ILM_kernel_variance::Array{Float64, 2},
      phylogenetic_kernel_variance::Array{Float64, 2},
      event_variance::Float64,
      event_extents::EventExtents,
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
               n::Int64,
               thin::Int64,
               ILM_kernel_variance::Array{Float64, 2},
               phylogenetic_kernel_variance::Array{Float64, 2},
               event_variance::Float64,
               event_extents::EventExtents,
               event_obs::EventObservations,
               seq_obs::Dict{Int64, Sequence},
               riskparameter_priors::RiskParameterPriors,
               riskfuncs::RiskFunctions,
               substitutionmodel_priors::SubstitutionModelPrior,
               population::DataFrame,
               iterprob::Vector{Float64},
               acceptance_rates=false::Bool)
  if thin < 0
    error("Thinning rate can not be less than 1")
  end
  progressbar = Progress(n, 5, "Performing $n iterations", 25)
  individuals = event_obs.individuals
  if typeof(event_obs) == SEIR_EventObservations
    validevents = find([.!isnan.(event_obs.infected) .!isnan.(event_obs.infected) .!isnan.(event_obs.removed)])
    eventdims = (individuals, 3)
  elseif typeof(event_obs) == SIR_EventObservations
    validevents = find([.!isnan.(event_obs.infected) .!isnan.(event_obs.removed)])
    eventdims = (individuals, 2)
  elseif typeof(event_obs) == SEI_EventObservations
    validevents = find([.!isnan.(event_obs.infected) .!isnan.(event_obs.infected)])
    eventdims = (individuals, 2)
  elseif typeof(event_obs) == SI_EventObservations
    validevents = find(.!isnan.(event_obs.infected))
    eventdims = (individuals, 1)
  end
  acceptance_rates_array = fill(0, (2, length(iterprob)))
  validexposures = find(.!isnan.(event_obs.infected))
  riskparameters_previous = copy(pathogen_trace.riskparameters[end])
  substitutionmodel_previous = copy(pathogen_trace.substitutionmodel[end])
  events_previous = copy(pathogen_trace.events[end])
  network_previous = copy(pathogen_trace.network[end])
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
      substitutionmodel_proposal = propose(pathogen_trace.substitutionmodel[end],
                                           phylogenetic_kernel_variance)
    else
      substitutionmodel_proposal = copy(substitutionmodel_previous)
    end

    if iterationtype == 3
      j, k = ind2sub(eventdims, sample(validevents))
      events_proposal = propose(j, k,
                                events_previous,
                                network_previous,
                                event_obs,
                                event_variance,
                                event_extents)
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
        j = sample(validexposures)
        network_proposal = propose(j,
                                   network_previous,
                                   network_rates)
      elseif iterationtype == 5
        network_proposal = propose(network_rates)
      else
        network_proposal = copy(network_previous)
      end

      tree_proposal, tree_ids, leaf_ids = generate_tree(events_proposal,
                                                        event_obs,
                                                        network_proposal)
      llikelihood += loglikelihood(tree_proposal,
                                   tree_ids,
                                   leaf_ids,
                                   substitutionmodel_proposal,
                                   seq_obs)
    else
      llikelihood = -Inf
    end

    lposterior = lprior + llikelihood

    if acceptance_rates
      acceptance_rates_array[1, iterationtype] += 1
    end
    if MHaccept(lposterior, lposterior_previous)
      riskparameters_previous = riskparameter_proposal
      substitutionmodel_previous = substitutionmodel_proposal
      events_previous = events_proposal
      network_previous = network_proposal
      lposterior_previous = lposterior
      if acceptance_rates
        acceptance_rates_array[2, iterationtype] += 1
      end
    end
    if mod(i, thin) == 0
      push!(pathogen_trace, PathogenIteration(copy(riskparameters_previous),
                                              copy(substitutionmodel_previous),
                                              copy(events_previous),
                                              copy(network_previous),
                                              copy(lposterior_previous)))
    end
  end
  if acceptance_rates
    return pathogen_trace, acceptance_rates_array
  else
    return pathogen_trace
  end
end
