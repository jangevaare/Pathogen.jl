# Partition without dropping the remainder when the event count is not divisible.
_event_batch(n::Integer, batches::Integer, batch::Integer) =
  (fld(n * (batch-1), batches)+1):fld(n * batch, batches)

function update!(rng::AbstractRNG, mc::MarkovChain{T, M},
                 mcmc::MCMC{T, M}, Σ::Array{Float64, 2}, σ::Float64;
                 condition_on_network::Bool=false, event_batches::Int64=1,
                 transmission_rate_cache::Union{Nothing, TransmissionRateCache}=nothing,
                 likelihood_workspace::Union{Nothing, LikelihoodWorkspace{T}}=nothing) where {
                 T <: DiseaseStateSequence, M <: TNILM}
  work = likelihood_workspace === nothing ? LikelihoodWorkspace(T, individuals(mcmc.population)) : likelihood_workspace
  new_terms = _ensure_network_probabilities!(mc, mcmc, work, transmission_rate_cache)
  new_tr_cache = transmission_rate_cache
  new_events = mc.events[end]
  new_rp = mc.risk_parameters[end]
  new_network = mc.transmission_network[end]
  new_base = mc.log_base_posterior
  event_indices = findall(convert(Vector{Float64}, new_events) .> -Inf)
  aug_order = sample(rng, event_indices, length(event_indices), replace=false)
  event_batches < 0 && throw(ArgumentError("Cannot have negative amount of event batches"))
  event_batches = min(event_batches, length(aug_order))
  rates_match = false
  local proposed_tr, proposed_tx
  for i = 1:(1 + event_batches)
    rates_match = false
    log_proposal_ratio = 0.0
    if i == 1
      proposed_events = new_events
      proposed_rp = generate(rng, RiskParameters{T}, new_rp, Σ)
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
    log_uniform = log(rand(rng))
    proposed_base = -Inf
    proposed_terms = new_terms
    log_acceptance = -Inf
    proposed_prior = logprior(proposed_rp, mcmc.risk_priors)
    if isfinite(proposed_prior) && isfinite(log_proposal_ratio)
      likelihood, proposed_tr, proposed_tx = loglikelihood(proposed_rp, mcmc.risk_functions,
        proposed_events, mcmc.population, mcmc.starting_states;
        transmission_rate_cache=proposed_tr_cache, workspace=work)
      if isfinite(likelihood) && _has_network_support(proposed_tr, mcmc.transmission_network_prior, proposed_tx)
        proposed_base = proposed_prior + likelihood
        log_acceptance = proposed_base - new_base + log_proposal_ratio
        if condition_on_network
          proposed_terms = _network_logprobabilities(new_network, proposed_tr,
            mcmc.transmission_network_prior, proposed_tx, mcmc.starting_states)
          log_acceptance += _retained_network_logratio(proposed_terms, new_terms, Int[])
        end
      end
    end
    if isfinite(proposed_base) && isfinite(log_acceptance) && log_uniform <= log_acceptance
      @debug "MCMC proposal accepted" log_acceptance
      rates_match = true
      new_rp, new_events, new_tr_cache = proposed_rp, proposed_events, proposed_tr_cache
      new_base, new_terms = proposed_base, proposed_terms
    else
      @debug "MCMC proposal rejected" log_acceptance
    end
  end
  if rates_match
    rates, transmissions = proposed_tr, proposed_tx
  else
    _, rates, transmissions = loglikelihood(new_rp, mcmc.risk_functions, new_events,
      mcmc.population, mcmc.starting_states; loglikelihood_output=false,
      transmission_rate_cache=new_tr_cache, workspace=work)
  end
  # A full conditional network draw is a Gibbs step in either TNILM mode.
  new_network = generate(rng, TransmissionNetwork, rates, mcmc, transmissions)
  new_terms = _network_logprobabilities(new_network, rates, mcmc.transmission_network_prior,
                                        transmissions, mcmc.starting_states)
  isfinite(sum(new_terms)) || error("Sampled network has zero conditional probability")
  _retain_cache!(transmission_rate_cache, new_tr_cache)
  mc.log_base_posterior = new_base
  mc.network_logprobabilities = new_terms
  push!(mc.events, new_events)
  push!(mc.transmission_network, new_network)
  push!(mc.risk_parameters, new_rp)
  push!(mc.log_posterior, new_base + sum(new_terms))
  mc.iterations += 1
  return mc
end

function update!(rng::AbstractRNG, mc::MarkovChain{T, M},
                 mcmc::MCMC{T, M}, Σrp::Array{Float64, 2}, Σsm::Array{Float64, 2}, σ::Float64;
                 event_batches::Int64=1,
                 transmission_rate_cache::Union{Nothing, TransmissionRateCache}=nothing,
                 DA_couple::Bool=false,
                 likelihood_workspace::Union{Nothing, LikelihoodWorkspace{T}}=nothing) where {
                 T <: DiseaseStateSequence, M <: PhyloILM}
  work = likelihood_workspace === nothing ? LikelihoodWorkspace(T, individuals(mcmc.population)) : likelihood_workspace
  new_terms = _ensure_network_probabilities!(mc, mcmc, work, transmission_rate_cache)
  new_tr_cache = transmission_rate_cache
  new_events, new_rp = mc.events[end], mc.risk_parameters[end]
  new_sm, new_network = mc.substitution_model[end], mc.transmission_network[end]
  new_base = mc.log_base_posterior
  event_indices = findall(convert(Vector{Float64}, new_events) .> -Inf)
  aug_order = sample(rng, event_indices, length(event_indices), replace=false)
  event_batches < 0 && throw(ArgumentError("Cannot have negative amount of event batches"))
  event_batches = min(event_batches, length(aug_order))
  steps = DA_couple ? 1 + event_batches : 2 + event_batches
  for i = 1:steps
    log_proposal_ratio = 0.0
    redrawn = Int[]
    time_batch = 1 < i < event_batches + 2
    if i == 1
      proposed_events = new_events
      proposed_rp = generate(rng, RiskParameters{T}, new_rp, Σrp)
      proposed_sm = generate(rng, mcmc.substitution_model, new_sm, Σsm)
      proposed_tr_cache = _proposal_cache!(work, individuals(new_events))
    elseif time_batch
      batch_indices = aug_order[_event_batch(length(aug_order), event_batches, i-1)]
      proposed_events, log_forward = _propose_events(rng, new_events, batch_indices,
        σ, mcmc.event_extents, mcmc.event_observations, new_network)
      if DA_couple
        redrawn = [index for index in batch_indices if index <= individuals(new_events)]
      end
      proposed_rp, proposed_sm, proposed_tr_cache = new_rp, new_sm, new_tr_cache
    else
      proposed_events, proposed_rp = new_events, new_rp
      proposed_sm, proposed_tr_cache = new_sm, new_tr_cache
    end
    log_uniform = log(rand(rng))
    proposed_base = -Inf
    proposed_terms = new_terms
    proposed_network = new_network
    log_acceptance = -Inf
    proposed_prior = logprior(proposed_rp, mcmc.risk_priors) +
                     logprior(proposed_sm, mcmc.substitution_model_prior)
    if isfinite(proposed_prior)
      likelihood, rates, transmissions = loglikelihood(proposed_rp, mcmc.risk_functions,
        proposed_events, mcmc.population, mcmc.starting_states;
        transmission_rate_cache=proposed_tr_cache, workspace=work)
      if isfinite(likelihood) && _has_network_support(rates, mcmc.transmission_network_prior, transmissions)
        if !DA_couple && i == steps && !isempty(transmissions)
          redrawn = [sample(rng, transmissions)]
        end
        if !isempty(redrawn)
          proposed_network = generate(rng, TransmissionNetwork, new_network, rates, mcmc, redrawn)
        end
        if time_batch
          log_proposal_ratio = _event_logproposal(new_events, proposed_events, batch_indices,
            σ, mcmc.event_extents, mcmc.event_observations, proposed_network) - log_forward
        end
        proposed_terms = _network_logprobabilities(proposed_network, rates,
          mcmc.transmission_network_prior, transmissions, mcmc.starting_states)
        if isfinite(sum(proposed_terms)) && isfinite(log_proposal_ratio)
          proposed_base = proposed_prior + likelihood
          # Source probabilities cancel only for parents actually redrawn. This
          # includes the reverse network proposal under the old event history.
          correction = _retained_network_logratio(proposed_terms, new_terms, redrawn) + log_proposal_ratio
          log_acceptance = proposed_base - new_base + correction
          # Sequence log probabilities are nonpositive, so this is an upper bound.
          if isfinite(log_acceptance) && log_uniform <= log_acceptance
            tree, obs_nodes = generate(Tree, proposed_events, mcmc.event_observations, proposed_network)
            leaf_data = Dict{Int64, RNASeq}()
            for k in eachindex(obs_nodes)
              if !isnothing(obs_nodes[k])
                leaf_data[obs_nodes[k]] = mcmc.event_observations.seq[k]
              end
            end
            proposed_base += PhyloModels.loglikelihood(tree, proposed_sm, leaf_data)
            log_acceptance = proposed_base - new_base + correction
          else
            proposed_base = -Inf
          end
        end
      end
    end
    if isfinite(proposed_base) && isfinite(log_acceptance) && log_uniform <= log_acceptance
      @debug "MCMC proposal accepted" log_acceptance
      new_rp, new_sm, new_events = proposed_rp, proposed_sm, proposed_events
      new_network, new_terms = proposed_network, proposed_terms
      new_base, new_tr_cache = proposed_base, proposed_tr_cache
    else
      @debug "MCMC proposal rejected" log_acceptance
    end
  end
  _retain_cache!(transmission_rate_cache, new_tr_cache)
  mc.log_base_posterior = new_base
  mc.network_logprobabilities = new_terms
  push!(mc.events, new_events)
  push!(mc.transmission_network, new_network)
  push!(mc.risk_parameters, new_rp)
  push!(mc.substitution_model, new_sm)
  push!(mc.log_posterior, new_base + sum(new_terms))
  mc.iterations += 1
  return mc
end

update!(mc::MarkovChain, args...; kwargs...) = update!(_default_rng(), mc, args...; kwargs...)
