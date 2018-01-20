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
               parameter_variance::Vector{Float64},
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
  riskparameter_params = length(riskparameter_priors)
  substitutionmodel_params = length(pathogen_trace.substitutionmodel[1].Θ)
  if n < 1
    error("The number of iterations must be > 0")
  elseif mod(n, thin) != 0
    error("The number of iterations must be a multiple of the thinning rate")
  elseif (riskparameter_params + substitutionmodel_params) != length(parameter_variance)
    error("The number of transition kernel variances provided does not match with the total number of model parameters")
  end

  progressbar = Progress(n, 5, "Performing $n iterations", 20)
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
  acceptance_rates_array = fill(0, (2, riskparameter_params + substitutionmodel_params + 2))
  validexposures = find(.!isnan.(event_obs.infected))
  riskparameters_previous = copy(pathogen_trace.riskparameters[end])
  substitutionmodel_previous = copy(pathogen_trace.substitutionmodel[end])
  events_previous = copy(pathogen_trace.events[end])
  network_previous = copy(pathogen_trace.network[end])
  lposterior_previous = copy(pathogen_trace.logposterior[end])

  # One substep for parameter updates
  # One substep for exposure network update
  # One substep for each eventtime to be augmented
  o = riskparameter_params + substitutionmodel_params + length(validevents) + 1

  for i = 1:n
    next!(progressbar, showvalues = [("Parameter proposals", "$(acceptance_rates_array[2, 1:(riskparameter_params + substitutionmodel_params)])/$(acceptance_rates_array[1, 1])");
                                     ("Event proposals", "$(acceptance_rates_array[2, end-1])/$(acceptance_rates_array[1, end-1])");
                                     ("Network proposals", "$(acceptance_rates_array[2, end])/$(acceptance_rates_array[1, end])")])
    parameter_order = sample(1:(riskparameter_params + substitutionmodel_params),
                             (riskparameter_params + substitutionmodel_params),
                             replace=false)
    augmentation_order = sample(validevents, length(validevents), replace=false)
    for j = 1:o
      riskparameter_proposal = copy(riskparameters_previous)
      substitutionmodel_proposal = copy(substitutionmodel_previous)
      if j <= riskparameter_params + substitutionmodel_params
        k = parameter_order[j]
        if k <= riskparameter_params
          riskparameter_proposal[k] = rand(Normal(riskparameter_proposal[k],
                                                  parameter_variance[k]))
        else
          substitutionmodel_proposal.Θ[k-riskparameter_params] = rand(Normal(substitutionmodel_proposal.Θ[k-riskparameter_params],
                                                                             parameter_variance[k]))
        end
      end

      if (riskparameter_params + substitutionmodel_params) < j & j < o
        k = augmentation_order[j - (riskparameter_params + substitutionmodel_params)]
        l, m = ind2sub(eventdims, k)
        events_proposal = propose(l, m,
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
        if j == o
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
            llikelihood += -Inf
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

      if j <= (riskparameter_params + substitutionmodel_params)
        acceptance_rates_array[1, parameter_order[j]] += 1
      elseif j < o
        acceptance_rates_array[1, end-1] += 1
      else
        acceptance_rates_array[1, end] += 1
      end

      if MHaccept(lposterior, lposterior_previous)
        riskparameters_previous = riskparameter_proposal
        substitutionmodel_previous = substitutionmodel_proposal
        events_previous = events_proposal
        network_previous = network_proposal
        lposterior_previous = lposterior
        if j <= (riskparameter_params + substitutionmodel_params)
          acceptance_rates_array[2, parameter_order[j]] += 1
        elseif j < o
          acceptance_rates_array[2, end-1] += 1
        else
          acceptance_rates_array[2, end] += 1
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
