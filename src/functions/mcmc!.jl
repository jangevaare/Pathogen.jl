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
      population::DataFrame;
      acceptance_rates=false::Bool,
      conditional_network_proposals=true::Bool)

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
               population::DataFrame;
               acceptance_rates=false::Bool,
               conditional_network_proposals=true::Bool)
  if n < 1
    error("The number of iterations must be > 0")
  elseif mod(n, thin) != 0
    error("The number of iterations must be a multiple of the thinning rate")
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
  acceptance_rates_array = fill(0, (2, 3))
  validexposures = find(.!isnan.(event_obs.infected))
  riskparameters_previous = copy(pathogen_trace.riskparameters[end])
  substitutionmodel_previous = copy(pathogen_trace.substitutionmodel[end])
  events_previous = copy(pathogen_trace.events[end])
  network_previous = copy(pathogen_trace.network[end])
  lposterior_previous = copy(pathogen_trace.logposterior[end])

  # One substep for parameter updates
  # One substep for exposure network update
  # One substep for each eventtime to be augmented
  m = length(validevents) + 1 + 1

  for i = 1:n
    next!(progressbar)
    augmentation_order = sample(validevents, length(validevents), replace=false)
    for j = 1:m
      if j == 1
        riskparameter_proposal = propose(riskparameters_previous,
                                         ILM_kernel_variance)
        substitutionmodel_proposal = propose(substitutionmodel_previous,
                                             phylogenetic_kernel_variance)
      else
        riskparameter_proposal = copy(riskparameters_previous)
        substitutionmodel_proposal = copy(substitutionmodel_previous)
      end

      if 1 < j & j < m
        k, l = ind2sub(eventdims, augmentation_order[j-1])
        events_proposal = propose(k, l,
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
        if j == m
          # Only sample individuals which have the possibility of having an internal
          # exoposure (assuming all have the possibility of external exposure)
          candidates = find(sum(network_rates.internal, 1) .> 0)
          if length(candidates) > 0
          # 95% of the time, propose a single change to tree
            if rand() < 0.95
              k = sample(candidates)
              network_proposal = propose(k,
                                         network_previous,
                                         network_rates,
                                         conditional_network_proposals = conditional_network_proposals)
            else
              network_proposal = propose(network_rates,
                                         conditional_network_proposals = conditional_network_proposals)
            end
          else
            network_proposal = copy(network_previous)
          end
        else
          network_proposal = copy(network_previous)
        end

        # Find loglikelihood of transmission network if rates from ILM were not
        # used in the generation of the proposal
        if !conditional_network_proposals
          llikelihood += loglikelihood(network_proposal,
                                       network_rates)
        end

        # Generate a tree
        tree_proposal = generate_tree(events_proposal,
                                      event_obs,
                                      network_proposal)

        # Calculate the tree's loglikelihood
        llikelihood += loglikelihood(tree_proposal,
                                     substitutionmodel_proposal,
                                     seq_obs)
      else
        llikelihood = -Inf
      end

      lposterior = lprior + llikelihood

      if acceptance_rates
        if j == 1
          acceptance_rates_array[1, 1] += 1
        elseif 1 < j & j < m
          acceptance_rates_array[1, 2] += 1
        elseif j == m
          acceptance_rates_array[1, 3] += 1
        end
      end
      if MHaccept(lposterior, lposterior_previous)
        riskparameters_previous = riskparameter_proposal
        substitutionmodel_previous = substitutionmodel_proposal
        events_previous = events_proposal
        network_previous = network_proposal
        lposterior_previous = lposterior
        if acceptance_rates
          if j == 1
            acceptance_rates_array[2, 1] += 1
          elseif 1 < j & j < m
            acceptance_rates_array[2, 2] += 1
          elseif j == m
            acceptance_rates_array[2, 3] += 1
          end
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
  end
  if acceptance_rates
    return pathogen_trace, acceptance_rates_array
  else
    return pathogen_trace
  end
end
