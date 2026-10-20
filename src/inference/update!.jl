# Partition without dropping the remainder when the event count is not divisible.
_event_batch(n::Integer, batches::Integer, batch::Integer) =
  (fld(n * (batch-1), batches)+1):fld(n * batch, batches)

function update!(rng::AbstractRNG, mc::MarkovChain{T, M},
                 mcmc::MCMC{T, M},
                 Σ::Array{Float64, 2},
                 σ::Float64;
                 condition_on_network::Bool=false,
                 event_batches::Int64=1,
                 transmission_rate_cache::Union{Nothing, TransmissionRateCache}=nothing,
                 likelihood_workspace::Union{Nothing, LikelihoodWorkspace{T}}=nothing) where {
                 T <: DiseaseStateSequence,
                 M <: TNILM}
  # Initialize
  work = likelihood_workspace === nothing ? LikelihoodWorkspace(T, individuals(mcmc.population)) : likelihood_workspace
  new_tr_cache = transmission_rate_cache
  new_events = mc.events[end]
  new_rp = mc.risk_parameters[end]
  new_network = mc.transmission_network[end]
  new_lposterior = mc.log_posterior[end]
  # Randomize event time augmentation
  event_indices = findall(convert(Vector{Float64}, new_events) .> -Inf)
  aug_order = sample(rng, event_indices, length(event_indices), replace=false)
  @debug "$(length(aug_order)) events to augment"
  if event_batches < 0
    throw(ArgumentError("Cannot have negative amount of event batches"))
  end
  if event_batches > length(aug_order)
    @warn "More event batches than there are events to augment ($(event_batches) > $(length(aug_order))), setting to maximum ($(length(aug_order)))"
    event_batches = length(aug_order)
  end
  @debug "Performing data augmentation in $event_batches batches"
  rates_match = false
  local proposed_tnr, proposed_tx
  for i = 1:(1 + event_batches)
    log_proposal_ratio = 0.0
    rates_match = false
    if i == 1
      # Propose new risk parameters
      proposed_events = new_events
      proposed_rp = generate(rng, RiskParameters{T}, new_rp, Σ)
      @debug "Generating new TransmissionRateCache for TN-ILM parameter proposal"
      proposed_tr_cache = _proposal_cache!(work, individuals(new_events))
    else
      batch_indices = aug_order[_event_batch(length(aug_order), event_batches, i-1)]
      network = condition_on_network ? new_network : nothing
      proposed_events, log_forward = _propose_events(rng, new_events, batch_indices,
        σ, mcmc.event_extents, mcmc.event_observations, network)
      log_proposal_ratio = _event_logproposal(new_events, proposed_events, batch_indices,
        σ, mcmc.event_extents, mcmc.event_observations, network) - log_forward
      proposed_rp = new_rp
      proposed_tr_cache = new_tr_cache
    end
    proposed_lprior = logprior(proposed_rp, mcmc.risk_priors)
    # Draw the Metropolis decision once; evaluate the complete epidemic likelihood.
    ll_acceptance_threshold = log(rand(rng)) + new_lposterior - proposed_lprior
    proposed_lposterior = -Inf
    if ll_acceptance_threshold < Inf && isfinite(log_proposal_ratio)
      proposed_llikelihood, proposed_tnr, proposed_tx = loglikelihood(proposed_rp,
                                           mcmc.risk_functions,
                                           proposed_events,
                                           mcmc.population,
                                           mcmc.starting_states,
                                           transmission_rate_cache = proposed_tr_cache, workspace = work)
      proposed_lposterior = proposed_lprior + proposed_llikelihood
    else
      proposed_llikelihood = -Inf
      proposed_lposterior = -Inf
    end
    if isfinite(proposed_llikelihood) && isfinite(log_proposal_ratio) && proposed_llikelihood + log_proposal_ratio >= ll_acceptance_threshold
      @debug "MCMC proposal accepted (acceptance probability = $(round(min(1.0, exp(proposed_lposterior - new_lposterior + log_proposal_ratio)), digits=3)))"
      rates_match = true
      new_rp = proposed_rp
      new_tr_cache = proposed_tr_cache
      new_events = proposed_events
      new_lposterior = proposed_lposterior
    else
      @debug "next!: MCMC proposal rejected (acceptance probability = $(round(min(1.0, exp(proposed_lposterior - new_lposterior + log_proposal_ratio)), digits=3)))"
    end
  end
  mc.iterations += 1
  if rates_match
    tnr, tx = proposed_tnr, proposed_tx
  else
    _, tnr, tx = loglikelihood(new_rp, mcmc.risk_functions, new_events,
                               mcmc.population, mcmc.starting_states;
                               loglikelihood_output=false,
                               transmission_rate_cache=new_tr_cache, workspace=work)
  end
  new_network = generate(rng, TransmissionNetwork, tnr, mcmc, tx)
  _retain_cache!(transmission_rate_cache, new_tr_cache)
  push!(mc.events, new_events)
  push!(mc.transmission_network, new_network)
  push!(mc.risk_parameters, new_rp)
  push!(mc.log_posterior, new_lposterior)
  return mc
end

function update!(rng::AbstractRNG, mc::MarkovChain{T, M},
                 mcmc::MCMC{T, M},
                 Σrp::Array{Float64, 2},
                 Σsm::Array{Float64, 2},
                 σ::Float64;
                 event_batches::Int64=1,
                 transmission_rate_cache::Union{Nothing, TransmissionRateCache}=nothing,
                 DA_couple::Bool=false,
                 likelihood_workspace::Union{Nothing, LikelihoodWorkspace{T}}=nothing) where {
                 T <: DiseaseStateSequence,
                 M <: PhyloILM}
  # Initialize
  work = likelihood_workspace === nothing ? LikelihoodWorkspace(T, individuals(mcmc.population)) : likelihood_workspace
  local tnr
  new_tr_cache = transmission_rate_cache
  new_events = mc.events[end]
  new_rp = mc.risk_parameters[end]
  new_sm = mc.substitution_model[end]
  new_network = mc.transmission_network[end]
  new_lposterior = mc.log_posterior[end]
  # Randomize event time augmentation
  event_indices = findall(convert(Vector{Float64}, new_events) .> -Inf)
  aug_order = sample(rng, event_indices, length(event_indices), replace=false)
  @debug "$(length(aug_order)) events to augment"
  if event_batches < 0
    throw(ArgumentError("Cannot have negative amount of event batches"))
  end
  if event_batches > length(aug_order)
    @warn "More event batches than there are events to augment ($(event_batches) > $(length(aug_order))), setting to maximum ($(length(aug_order)))"
    event_batches = length(aug_order)
  end
  @debug "Performing data augmentation in $event_batches batches"
  steps = DA_couple ? 1 + event_batches : 1 + event_batches + 1
  for i = 1:steps
    log_proposal_ratio = 0.0
    if i == 1
      # Propose new risk parameters
      proposed_events = new_events
      proposed_rp = generate(rng, RiskParameters{T}, new_rp, Σrp)
      proposed_sm = generate(rng, mcmc.substitution_model, new_sm, Σsm)
      @debug "Generating new TransmissionRateCache for TN-ILM parameter proposal"
      proposed_tr_cache = _proposal_cache!(work, individuals(new_events))
    elseif i < (event_batches + 2)
      batch_indices = aug_order[_event_batch(length(aug_order), event_batches, i-1)]
      proposed_events, log_forward = _propose_events(rng, new_events, batch_indices,
        σ, mcmc.event_extents, mcmc.event_observations, new_network)
      n_ids = individuals(new_events)
      t_ids = [index for index in batch_indices if index <= n_ids]
      proposed_rp = new_rp
      proposed_sm = new_sm
      proposed_tr_cache = new_tr_cache
    elseif i == (event_batches + 2)
      proposed_events = new_events
      proposed_rp = new_rp
      proposed_sm = new_sm
      proposed_tr_cache = new_tr_cache
    end
    proposed_lprior = logprior(proposed_rp, mcmc.risk_priors)
    proposed_lprior += logprior(proposed_sm, mcmc.substitution_model_prior)
    # Draw the Metropolis decision once; evaluate the complete epidemic likelihood.
    ll_acceptance_threshold = log(rand(rng)) + new_lposterior - proposed_lprior
    proposed_lposterior = -Inf
    if ll_acceptance_threshold < Inf && isfinite(log_proposal_ratio)
      proposed_llikelihood, tnr, tx = loglikelihood(proposed_rp,
                                                    mcmc.risk_functions,
                                                    proposed_events,
                                                    mcmc.population,
                                                    mcmc.starting_states,
                                                    transmissions_output = true,
                                                    transmission_rate_cache = proposed_tr_cache, workspace = work)
      if isfinite(proposed_llikelihood)
        if DA_couple
          if i == 1 || length(t_ids) == 0
            proposed_network = new_network
          elseif 1 < i < event_batches + 2
            proposed_network = generate(rng, TransmissionNetwork, new_network, tnr, mcmc, t_ids)
            @debug "Generating a Transmission Network proposal" ProposedNetwork = proposed_network
          end
        else
          if i < event_batches + 2
            proposed_network = new_network
          elseif i == event_batches + 2
            proposed_network = isempty(tx) ? new_network : generate(rng, TransmissionNetwork, new_network, tnr, mcmc, sample(rng, tx))
            @debug "Generating a transmission network proposal" ProposedNetwork = proposed_network
          end
        end
        if 1 < i < event_batches + 2
          log_proposal_ratio = _event_logproposal(new_events, proposed_events, batch_indices,
            σ, mcmc.event_extents, mcmc.event_observations, proposed_network) - log_forward
        end
        if isfinite(proposed_llikelihood) && isfinite(log_proposal_ratio) && proposed_llikelihood + log_proposal_ratio >= ll_acceptance_threshold
          @debug "Generating a Tree for proposed network and event times"
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
      end
    else
      proposed_llikelihood = -Inf
      proposed_lposterior = -Inf
    end
    if isfinite(proposed_llikelihood) && isfinite(log_proposal_ratio) && proposed_llikelihood + log_proposal_ratio >= ll_acceptance_threshold
      @debug "MCMC proposal accepted (acceptance probability = $(round(min(1.0, exp(proposed_lposterior - new_lposterior + log_proposal_ratio)), digits=3)))"
      new_rp = proposed_rp
      new_sm = proposed_sm
      new_tr_cache = proposed_tr_cache
      new_events = proposed_events
      new_network = proposed_network
      new_lposterior = proposed_lposterior
    else
      @debug "MCMC proposal rejected (acceptance probability = $(round(min(1.0, exp(proposed_lposterior - new_lposterior + log_proposal_ratio)), digits=3)))"
    end
  end
  _retain_cache!(transmission_rate_cache, new_tr_cache)
  mc.iterations += 1
  push!(mc.events, new_events)
  push!(mc.transmission_network, new_network)
  push!(mc.risk_parameters, new_rp)
  push!(mc.substitution_model, new_sm)
  push!(mc.log_posterior, new_lposterior)
  return mc
end

update!(mc::MarkovChain, args...; kwargs...) = update!(_default_rng(), mc, args...; kwargs...)
