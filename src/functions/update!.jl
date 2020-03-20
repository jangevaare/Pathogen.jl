function update!(events::Events{T},
                 event::Event{T}) where T <: EpidemicModel
  events[event.new_state][event.individual] = event.time
  return events
end

function update!(states::Vector{DiseaseState},
                 event::Event{T}) where T <: EpidemicModel
  states[event.individual] = event.new_state
  return states
end

function update!(net::TransmissionNetwork,
                 xfer::ExogenousTransmission)
  net.external[xfer.individual] = true
  return net
end

function update!(net::TransmissionNetwork,
                 xfer::EndogenousTransmission)
  net.internal[xfer.source, xfer.individual] = true
  return net
end

function update!(net::TransmissionNetwork,
                 xfer::NoTransmission)
  return net
end

function update!(tr::TransmissionRates,
                 event::Event{T},
                 states::Vector{DiseaseState},
                 pop::Population,
                 rf::RiskFunctions{T},
                 rp::RiskParameters{T};
                 xzero::Bool=false) where T <: EpidemicModel
  id = event.individual
  if _new_transmission(event) && xzero
    tr.external[id] = 0.0
    tr.internal[:, id] .= 0.0
  end
  if event.new_state == State_I
    transmissibility = rf.transmissibility(rp.transmissibility, pop, id)
    @simd for i in findall(states .== Ref(State_S))
      tr.internal[id, i] = rf.susceptibility(rp.susceptibility, pop, i) *
                           rf.infectivity(rp.infectivity, pop, i, id) *
                           transmissibility
    end
  elseif event.new_state == State_R
    @simd for i in findall(states .== Ref(State_S))
      tr.internal[id, i] = 0.0
    end
  end
  return tr
end

function update!(rates::EventRates{T},
                 tr::TransmissionRates,
                 event::Event{T},
                 states::Vector{DiseaseState},
                 pop::Population,
                 rf::RiskFunctions{T},
                 rp::RiskParameters{T}) where T <: EpidemicModel
  id = event.individual
  if event.new_state == State_E
    rates.exposure[id] = 0.0
    rates.infection[id] = rf.latency(rp.latency, pop, id)
    @logmsg LogLevel(-5000) "Exposure rate total for i = $id updated"  λ = rates.exposure[id]
    @logmsg LogLevel(-5000) "Infection rate for i = $id updated" λ = rates.infection[id]
  elseif event.new_state == State_I
    rates.infection[id] = 0.0
    @logmsg LogLevel(-5000) "Infection rate for i = $id updated" λ = rates.infection[id]
    if T in [SEIR; SIR]
      rates.removal[id] = rf.removal(rp.removal, pop, id)
      @logmsg LogLevel(-5000) "Removal rate for i = $id updated" λ = rates.removal[id]
    end
    if T in [SEIR; SEI]
      @simd for i in findall(states .== Ref(State_S))
        # This assumes `TransmissionRates` already updated!
        rates.exposure[i] += tr.internal[id, i]
        @logmsg LogLevel(-5000) "Exposure rate total for i = $id updated"  λ = rates.exposure[id]
      end
    elseif T in [SIR; SI]
      @simd for i in findall(states .== Ref(State_S))
        # This assumes `TransmissionRates` already updated!
        rates.infection[i] += tr.internal[id, i]
        @logmsg LogLevel(-5000) "Infection rate total for i = $id updated"  λ = rates.infection[id]
      end
    end
  elseif event.new_state == State_R
    rates.removal[id] = 0.0
    @logmsg LogLevel(-5000) "Removal rate for i = $id updated" λ = rates.removal[id]
    if T == SEIR
      @simd for i in findall(states .== Ref(State_S))
        # This assumes `TransmissionRates` already updated!
        rates.exposure[i] = sum(tr.internal[:, i])
        @logmsg LogLevel(-5000) "Exposure rate total for i = $i updated"  λ = rates.exposure[i]
      end
    elseif T == SIR
      @simd for i in findall(states .== Ref(State_S))
        # This assumes `TransmissionRates` already updated!
        rates.infection[i] = sum(tr.internal[:, i])
        @logmsg LogLevel(-5000) "Infection rate for i = $i updated" λ = rates.infection[i]
      end
    end
  end
  return rates
end

function update!(tr::TransmissionRates,
                 rates::EventRates{T},
                 event::Event{T},
                 states::Vector{DiseaseState},
                 pop::Population,
                 rf::RiskFunctions{T},
                 rp::RiskParameters{T}) where T <: EpidemicModel
  update!(tr, event, states, pop, rf, rp)
  update!(rates, tr, event, states, pop, rf, rp)
  return tr, rates
end

function update!(s::Simulation{T}, event::AbstractEvent) where T <: EpidemicModel
  if event.time < s.time
    @error "Time of a new event must be >= that of the previous event"
  elseif event.time < Inf
    update!(s.events, event)
    update!(s.disease_states, event)
    update!(s.transmission_network,
            generate(Transmission,
                     s.transmission_rates,
                     event))
    update!(s.transmission_rates,
            s.event_rates,
            event,
            s.disease_states,
            s.population,
            s.risk_functions,
            s.risk_parameters)
    # Count the iteration as having occurred.
    s.iterations += 1
  end
  # Update simulation time
  s.time = event.time
  # Return updated Simulation object
  return s
end

function update!(mc::MarkovChain{T},
                 mcmc::MCMC{T},
                 Σ::Array{Float64, 2},
                 σ::Float64;
                 condition_on_network::Bool=false,
                 event_batches::Int64=1) where T <: EpidemicModel
  # Initialize
  new_events = mc.events[end]
  new_events_array = new_events[_state_progressions[T][2:end]]
  new_params = mc.risk_parameters[end]
  new_network = mc.transmission_network[end]
  new_lposterior = mc.log_posterior[end]
  # Randomize event time augmentation
  event_indices = findall(new_events_array[:] .> -Inf)
  aug_order = sample(event_indices, length(event_indices), replace=false)
  @debug "$(length(aug_order)) events to augment"
  if event_batches < 0
    @error "Cannot have negative amount of event batches"
  end
  if event_batches > length(aug_order)
    @warn "More event batches than there are events to augment ($(event_batches) > $(length(aug_order))), setting to maximum ($(length(aug_order)))"
    event_batches = length(aug_order)
  end
  batch_size = fld(length(aug_order), event_batches)
  @debug "Performing data augmentation in batches of $batch_size events at a time"
  for i = 1:(event_batches + 1)
    if i <= event_batches
      for j = (batch_size*(i-1) + 1):minimum([(batch_size*i + 1); length(aug_order)])
        id, state_index = Tuple(CartesianIndices((new_events.individuals,
                                            length(_state_progressions[T][2:end])))[aug_order[j]])
        new_state = _state_progressions[T][state_index+1]
        time = new_events[new_state][id]
        # Conditioning on network means that only event times valid under current network will be proposed. This is useful for models which may have additional contributions to the posterior based on network, and require more modest proposals (e.g. phylodynamic models).
        if condition_on_network
          proposed_event = generate(Event,
                                    Event{T}(time, id, new_state),
                                    σ,
                                    mcmc.event_extents,
                                    mcmc.event_observations,
                                    new_events,
                                    new_network)
        else
          proposed_event = generate(Event,
                                    Event{T}(time, id, new_state),
                                    σ,
                                    mcmc.event_extents,
                                    mcmc.event_observations,
                                    new_events)
        end
        # For the first event of batch, we'll create an events array from current Markov chain position
        # For other events in batch, we'll update proposal itself
        if j == (batch_size*(i-1) + 1)
          proposed_events_array = reshape([new_events_array[1:(aug_order[j]-1)]
                                           proposed_event.time
                                           new_events_array[(aug_order[j]+1):end]],
                                          size(new_events_array))
        else
          proposed_events_array = reshape([proposed_events_array[1:(aug_order[j]-1)]
                                           proposed_event.time
                                           proposed_events_array[(aug_order[j]+1):end]],
                                          size(new_events_array))
        end
        # Only need to generate new `Events` on last event of batch
        if j == minimum([(batch_size*i + 1); length(aug_order)])
          proposed_events = Events{T}(proposed_events_array)
        end
      end
      proposed_params = new_params
    else
      # Propose new risk parameters
      proposed_events = new_events
      proposed_events_array = new_events_array
      proposed_params = generate(RiskParameters{T}, new_params, Σ)
    end
    proposed_lprior = logpriors(proposed_params, mcmc.risk_priors)
    # Based on the logprior and competiting MCMC iteration, this loglikelihood is required for acceptance
    # Calculating this in advance allows us to cut loglikelihood calculation short if it goes below threshold
    ll_acceptance_threshold = log(rand()) + new_lposterior - proposed_lprior
    if ll_acceptance_threshold < Inf
      proposed_llikelihood = loglikelihood(proposed_params,
                                           mcmc.risk_functions,
                                           proposed_events,
                                           mcmc.population,
                                           mcmc.starting_states,
                                           transmission_network_output = false,
                                           early_decision_value = ll_acceptance_threshold)
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
  new_network = loglikelihood(new_params,
                              mcmc.risk_functions,
                              new_events,
                              mcmc.population,
                              mcmc.starting_states,
                              loglikelihood_output = false,
                              transmission_network_output = true)
  push!(mc.events, new_events)
  push!(mc.transmission_network, new_network)
  push!(mc.risk_parameters, new_params)
  push!(mc.log_posterior, new_lposterior)
  return mc
end

function update!(mcmc::MCMC{T},
                 Σ::Array{Float64, 2},
                 σ::Float64) where T <: EpidemicModel
  @simd for mc in mcmc.markov_chains
    next!(mc, mcmc, Σ, σ)
  end
  return mcmc
end

