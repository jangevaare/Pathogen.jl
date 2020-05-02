function update!(mc::MarkovChain{S, M},
                 mcmc::MCMC{S, M},
                 Σ::Array{Float64, 2},
                 σ::Float64;
                 condition_on_network::Bool=false,
                 event_batches::Int64=1) where {
                 S <: DiseaseStateSequence, 
                 M <: TNILM}
  # Initialize
  new_events = mc.events[end]
  new_events_array = new_events[convert(DiseaseStates, S)[2:end]]
  new_params = mc.risk_parameters[end]
  new_network = mc.transmission_network[end]
  new_lposterior = mc.log_posterior[end]
  # Randomize event time augmentation
  event_indices = findall(new_events_array[:] .> -Inf)
  aug_order = sample(event_indices, length(event_indices), replace=false)
  @debug "$(length(aug_order)) events to augment"
  if event_batches < 0
    throw(ErrorException("Cannot have negative amount of event batches"))
  end
  if event_batches > length(aug_order)
    @warn "More event batches than there are events to augment ($(event_batches) > $(length(aug_order))), setting to maximum ($(length(aug_order)))"
    event_batches = length(aug_order)
  end
  batch_size = fld(length(aug_order), event_batches)
  @debug "Performing data augmentation in batches of $batch_size events at a time"
  for i = 1:(event_batches + 1)
    if i <= event_batches
      proposed_events_array = copy(new_events_array)
      proposed_events = Events{S}(proposed_events_array)
      for j = (batch_size*(i-1) + 1):min(batch_size*i + 1, length(aug_order))
        id, state_index = Tuple(CartesianIndices((individuals(new_events),
                                            length(convert(DiseaseStates, S)[2:end])))[aug_order[j]])
        new_state = convert(DiseaseStates, S)[state_index+1]
        time = new_events[new_state][id]
        # Conditioning on network means that only event times valid under current network will be proposed.
        # This is useful for models which may have additional contributions to the posterior based on network,
        # and require more modest proposals (e.g. phylodynamic models).
        if condition_on_network
          proposed_event = generate(AbstractEvent,
                                    Event{S}(time, id, new_state),
                                    σ,
                                    mcmc.event_extents,
                                    mcmc.event_observations,
                                    proposed_events,
                                    new_network)
        else
          proposed_event = generate(AbstractEvent,
                                    Event{S}(time, id, new_state),
                                    σ,
                                    mcmc.event_extents,
                                    mcmc.event_observations,
                                    proposed_events)
        end
        proposed_events_array = reshape([proposed_events_array[1:(aug_order[j]-1)]
                                         _time(proposed_event)
                                         proposed_events_array[(aug_order[j]+1):end]],
                                        size(new_events_array))
        proposed_events = Events{S}(proposed_events_array)
      end
      proposed_params = new_params
    else
      # Propose new risk parameters
      proposed_events = new_events
      proposed_events_array = new_events_array
      proposed_params = generate(RiskParameters{S}, new_params, Σ)
    end
    proposed_lprior = logprior(proposed_params, mcmc.risk_priors)
    # Based on the logprior and competiting MCMC iteration, this loglikelihood is required for acceptance
    # Calculating this in advance allows us to cut loglikelihood calculation short if it goes below threshold
    ll_acceptance_threshold = log(rand()) + new_lposterior - proposed_lprior
    if ll_acceptance_threshold < Inf
      proposed_llikelihood = loglikelihood(proposed_params,
                                           mcmc.risk_functions,
                                           proposed_events,
                                           mcmc.population,
                                           mcmc.starting_states,
                                           transmission_rates_output = false,
                                           transmissions_output = false,
                                           early_decision_value = ll_acceptance_threshold)[1]
      proposed_lposterior = proposed_lprior + proposed_llikelihood
    else
      proposed_llikelihood = -Inf
      proposed_lposterior = -Inf
    end
    if proposed_llikelihood >= ll_acceptance_threshold
      @debug "MCMC proposal accepted (acceptance probability = $(round(min(1.0, exp(proposed_lposterior - new_lposterior)), digits=3)))"
      new_params = proposed_params
      new_events = proposed_events
      new_events_array = proposed_events_array
      new_lposterior = proposed_lposterior
    else
      @debug "next!: MCMC proposal rejected (acceptance probability = $(round(min(1.0, exp(proposed_lposterior - new_lposterior)), digits=3)))"
    end
  end
  mc.iterations += 1
  tnr, tx = loglikelihood(new_params,
                          mcmc.risk_functions,
                          new_events,
                          mcmc.population,
                          mcmc.starting_states,
                          loglikelihood_output = false)[[2; 3]]
  new_network = generate(TransmissionNetwork, tnr, mcmc, tx)
  push!(mc.events, new_events)
  push!(mc.transmission_network, new_network)
  push!(mc.risk_parameters, new_params)
  push!(mc.log_posterior, new_lposterior)
  return mc
end

function update!(mcmc::MCMC{S, M},
                 Σ::Array{Float64, 2},
                 σ::Float64;
                 condition_on_network::Bool=false,
                 event_batches::Int64=1) where {
                 S <: DiseaseStateSequence,
                 M <: TNILM}
  @simd for mc in mcmc.markov_chains
    next!(mc, mcmc, Σ, σ, condition_on_network = condition_on_network, event_batches = event_batches)
  end
  return mcmc
end

function update!(mc::MarkovChain{S, M},
                 mcmc::MCMC{S, M},
                 Σrp::Array{Float64, 2},
                 Σsm::Array{Float64, 2},
                 σ::Float64;
                 condition_on_network::Bool=true,
                 event_batches::Int64=1) where {
                 S <: DiseaseStateSequence, 
                 M <: PhyloILM}
  # Initialize
  new_events = mc.events[end]
  new_events_array = new_events[convert(DiseaseStates, S)[2:end]]
  new_rp = mc.risk_parameters[end]
  new_sm = mc.substitution_model[end]
  new_network = mc.transmission_network[end]
  new_lposterior = mc.log_posterior[end]
  # Randomize event time augmentation
  event_indices = findall(new_events_array[:] .> -Inf)
  aug_order = sample(event_indices, length(event_indices), replace=false)
  @debug "$(length(aug_order)) events to augment"
  if event_batches < 0
    throw(ErrorException("Cannot have negative amount of event batches"))
  end
  if event_batches > length(aug_order)
    @warn "More event batches than there are events to augment ($(event_batches) > $(length(aug_order))), setting to maximum ($(length(aug_order)))"
    event_batches = length(aug_order)
  end
  batch_size = fld(length(aug_order), event_batches)
  @debug "Performing data augmentation in batches of $batch_size events at a time"
  for i = 1:(event_batches + 2)
    @debug "Beginning sub-step $i"
    if i <= event_batches
      for j = (batch_size*(i-1) + 1):minimum([(batch_size*i + 1); length(aug_order)])
        id, state_index = Tuple(CartesianIndices((individuals(new_events),
                                            length(convert(DiseaseStates, S)[2:end])))[aug_order[j]])
        new_state = convert(DiseaseStates, S)[state_index+1]
        time = new_events[new_state][id]
        if condition_on_network
          proposed_event = generate(AbstractEvent,
                                    Event{S}(time, id, new_state),
                                    σ,
                                    mcmc.event_extents,
                                    mcmc.event_observations,
                                    new_events,
                                    new_network)
        else
          proposed_event = generate(AbstractEvent,
                                    Event{S}(time, id, new_state),
                                    σ,
                                    mcmc.event_extents,
                                    mcmc.event_observations,
                                    new_events)
        end


        # For the first event of batch, we'll create an events array from current Markov chain position
        # For other events in batch, we'll update proposal itself
        if j == (batch_size*(i-1) + 1)
          proposed_events_array = reshape([new_events_array[1:(aug_order[j]-1)]
                                           _time(proposed_event)
                                           new_events_array[(aug_order[j]+1):end]],
                                          size(new_events_array))
        else
          proposed_events_array = reshape([proposed_events_array[1:(aug_order[j]-1)]
                                           _time(proposed_event)
                                           proposed_events_array[(aug_order[j]+1):end]],
                                          size(new_events_array))
        end
        # Only need to generate new `Events` on last event of batch
        if j == minimum([(batch_size*i + 1); length(aug_order)])
          proposed_events = Events{S}(proposed_events_array)
        end
      end
      if !condition_on_network
        proposed_network = nothing # generate later following ll calculation
      else
        proposed_network = new_network
      end
      proposed_rp = new_rp
      proposed_sm = new_sm
    elseif i == event_batches + 1 # Parameter update
      # Propose new risk parameters
      proposed_events = new_events
      proposed_events_array = new_events_array
      proposed_network = new_network
      proposed_rp = generate(RiskParameters{S}, new_rp, Σrp)
      proposed_sm = generate(mcmc.substitution_model, new_sm, Σsm)
    elseif i == event_batches + 2 # Network update
      proposed_events = new_events
      proposed_events_array = new_events_array
      proposed_network = nothing # generate later following ll calculation
      proposed_rp = new_rp
      proposed_sm = new_sm
    end
    proposed_lprior = logprior(proposed_rp, mcmc.risk_priors)
    proposed_lprior += logprior(proposed_sm, mcmc.substitution_model_prior)

    # Based on the logprior and competiting MCMC iteration, this loglikelihood is required for acceptance
    # Calculating this in advance allows us to cut loglikelihood calculation short if it goes below threshold
    ll_acceptance_threshold = log(rand()) + new_lposterior - proposed_lprior
    if ll_acceptance_threshold < Inf
      proposed_llikelihood, tnr, tx  = loglikelihood(proposed_rp,
                                                     mcmc.risk_functions,
                                                     proposed_events,
                                                     mcmc.population,
                                                     mcmc.starting_states,
                                                     transmission_rates_output = true,
                                                     transmissions_output = true,
                                                     early_decision_value = ll_acceptance_threshold)
      if proposed_llikelihood > ll_acceptance_threshold
        if (i == event_batches + 2) || (i <= event_batches&& !condition_on_network)
          proposed_network = generate(TransmissionNetwork, tnr, mcmc, tx)
          @debug "Generating a Transmission Network proposal" ProposedNetwork = proposed_network
        end
        @debug "Generating a Tree for proposed network and event times" ProposedNetwork = proposed_network
        tree, obs_nodes = generate(Tree, proposed_events, mcmc.event_observations, proposed_network)
        leaf_data = Dict{Int64, RNASeq}()
        for k = eachindex(obs_nodes)
          if !isnothing(obs_nodes[k])
            leaf_data[obs_nodes[k]] = mcmc.event_observations.seq[k]
          end
        end
        proposed_llikelihood += PhyloModels.loglikelihood(tree, proposed_sm, leaf_data)
        proposed_lposterior = proposed_lprior + proposed_llikelihood
      end
    else
      proposed_llikelihood = -Inf
      proposed_lposterior = -Inf
    end
    if proposed_llikelihood >= ll_acceptance_threshold
      @debug "MCMC proposal accepted (acceptance probability = $(round(min(1.0, exp(proposed_lposterior - new_lposterior)), digits=3)))"
      new_rp = proposed_rp
      new_sm = proposed_sm
      new_events = proposed_events
      new_events_array = proposed_events_array
      new_network = proposed_network
      new_lposterior = proposed_lposterior
    else
      @debug "next!: MCMC proposal rejected (acceptance probability = $(round(min(1.0, exp(proposed_lposterior - new_lposterior)), digits=3)))"
    end
  end
  mc.iterations += 1
  push!(mc.events, new_events)
  push!(mc.transmission_network, new_network)
  push!(mc.risk_parameters, new_rp)
  push!(mc.substitution_model, new_sm)
  push!(mc.log_posterior, new_lposterior)
  return mc
end

function update!(mcmc::MCMC{S, M},
                 Σrp::Array{Float64, 2},
                 Σsm::Array{Float64, 2},
                 σ::Float64;
                 condition_on_network::Bool=true,
                 event_batches::Int64=1) where {
                 S <: DiseaseStateSequence,
                 M <: PhyloILM}
  @simd for mc in mcmc.markov_chains
    next!(mc, mcmc, Σrp, Σsm, σ, condition_on_network = condition_on_network, event_batches = event_batches)
  end
  return mcmc
end
