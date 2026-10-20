function generate(rng::AbstractRNG, ::Type{Transmission}, tr::TransmissionRates,
                  prior::Union{Nothing, TNDistribution}, id::Int64)
  internal = [_transmission_weight(tr, prior, source, id) for source in 1:individuals(tr)]
  external_or_internal = Weights([_transmission_weight(tr, prior, 0, id), sum(internal)])
  if sum(external_or_internal) == 0.0
    @error "All transmission rates = 0.0, No transmission generated" external_or_internal
    return NoTransmission()
  elseif sample(rng, [true; false], external_or_internal)
    return ExogenousTransmission(id)
  else
    source = sample(rng, 1:individuals(tr), Weights(internal))
    return EndogenousTransmission(id, source)
  end
end

function generate(
  rng::AbstractRNG,
  ::Type{TransmissionNetwork},
  starting_states::DiseaseStates,
  tr::TransmissionRates,
  tnd::Union{Nothing, TNDistribution},
  ids::Vector{Int64})
  tn = TransmissionNetwork(starting_states)
  for i in ids
    tx = generate(rng, Transmission, tr, tnd, i)
    update!(tn, tx)
  end
  dropzeros!(tn.internal)
  return tn
end

function generate(
  rng::AbstractRNG,
  ::Type{TransmissionNetwork},
  tr::TransmissionRates,
  mcmc::MCMC,
  ids::Vector{Int64})
  return generate(rng, TransmissionNetwork, mcmc.starting_states, tr, mcmc.transmission_network_prior, ids)
end

function generate(
  rng::AbstractRNG,
  ::Type{TransmissionNetwork},
  lasttn::TransmissionNetwork,
  tr::TransmissionRates,
  tnd::Union{Nothing, TNDistribution},
  ids)
  tn = copy(lasttn)
  for i in ids
    tn.external[i] = 0
    for k in nzrange(tn.internal, i)
      nonzeros(tn.internal)[k] = false
    end
    tx = generate(rng, Transmission, tr, tnd, i)
    update!(tn, tx)
  end
  dropzeros!(tn.internal)
  return tn
end

function generate(
  rng::AbstractRNG,
  ::Type{TransmissionNetwork},
  lasttn::TransmissionNetwork,
  tr::TransmissionRates,
  mcmc::MCMC,
  ids)
  return generate(rng, TransmissionNetwork, lasttn, tr, mcmc.transmission_network_prior, ids)
end

function generate(
  rng::AbstractRNG,
  ::Type{Events},
  obs::EventObservations{T, M},
  extents::EventExtents{T}) where {
  T <: DiseaseStateSequence,
  M <: ILM}
  events = Events{T}(individuals(obs))
  exposed_state = State_E ∈ T
  removed_state = State_R ∈ T
  for i = 1:individuals(obs)
    if obs.infection[i] == -Inf
      update!(events, Event{T}(-Inf, i, State_I))
      if exposed_state
        update!(events, Event{T}(-Inf, i, State_E))
      end
      if removed_state && obs.removal[i] == -Inf
        update!(events, Event{T}(-Inf, i, State_R))
      end
    elseif !isnan(obs.infection[i])
      i_lb = obs.infection[i] - extents.infection[2]
      i_ub = obs.infection[i] - extents.infection[1]
      i_time = rand(rng, Uniform(i_lb, i_ub))
      update!(events, Event{T}(i_time, i, State_I))
      if exposed_state
        e_lb = i_time - extents.exposure[2]
        e_ub = i_time - extents.exposure[1]
        e_time = rand(rng, Uniform(e_lb, e_ub))
        update!(events, Event{T}(e_time, i, State_E))
      end
      if removed_state && !isnan(obs.removal[i])
        r_lb = maximum([obs.infection[i]; obs.removal[i] - extents.removal[2]])
        r_ub = obs.removal[i] - extents.removal[1]
        r_time = rand(rng, Uniform(r_lb, r_ub))
        update!(events, Event{T}(r_time, i, State_R))
      end
    end
  end
  return events
end

function generate(
  rng::AbstractRNG,
  ::Type{Events},
  mcmc::MCMC{T}) where {
  T <: DiseaseStateSequence}
  return generate(rng, Events, mcmc.event_observations, mcmc.event_extents)
end

function _event_bounds(event, extents, obs, events, network)
  return network === nothing ? _bounds(event, extents, obs, events) :
                               _bounds(event, extents, obs, events, network)
end

function _event_distribution(event, σ, extents, obs, events, network=nothing)
  lower, upper = _event_bounds(event, extents, obs, events, network)
  return truncated(Normal(_time(event), σ), lower, upper)
end

function generate(rng::AbstractRNG, ::Type{Event}, last_event::Event{T},
                  σ::Float64, extents::EventExtents{T}, obs::EventObservations,
                  events::Events{T}) where T <: DiseaseStateSequence
  return Event(rand(rng, _event_distribution(last_event, σ, extents, obs, events)), last_event)
end

function generate(rng::AbstractRNG, ::Type{Event}, last_event::Event{T},
                  σ::Float64, extents::EventExtents{T}, obs::EventObservations,
                  events::Events{T}, network::TransmissionNetwork) where T <: DiseaseStateSequence
  return Event(rand(rng, _event_distribution(last_event, σ, extents, obs, events, network)), last_event)
end

# A batch is a sequence of conditional proposals. Its reverse density uses the
# same index order, starting at the proposed history and restoring the old times.
function _propose_events(rng::AbstractRNG, current::Events{T}, indices,
                         σ::Float64, extents, obs, network=nothing) where T
  proposed = copy(current)
  σ == 0.0 && return proposed, 0.0
  states = convert(DiseaseStates, T)[2:end]
  n = individuals(current)
  log_forward = 0.0
  for index in indices
    id, state = mod1(index, n), states[fld(index-1, n)+1]
    event = Event{T}(proposed[state][id], id, state)
    distribution = _event_distribution(event, σ, extents, obs, proposed, network)
    time = rand(rng, distribution)
    log_forward += logpdf(distribution, time)
    update!(proposed, Event(time, event))
  end
  return proposed, log_forward
end

function _event_logproposal(target::Events{T}, source::Events{T}, indices,
                            σ::Float64, extents, obs, network=nothing) where T
  intermediate = copy(source)
  states = convert(DiseaseStates, T)[2:end]
  n = individuals(source)
  log_density = 0.0
  for index in indices
    id, state = mod1(index, n), states[fld(index-1, n)+1]
    event = Event{T}(intermediate[state][id], id, state)
    time = target[state][id]
    if σ == 0.0
      time == event.time || return -Inf
      continue
    end
    lower, upper = _event_bounds(event, extents, obs, intermediate, network)
    lower < upper && lower <= time <= upper || return -Inf
    log_density += logpdf(truncated(Normal(event.time, σ), lower, upper), time)
    update!(intermediate, Event(time, event))
  end
  return log_density
end

function generate(
  rng::AbstractRNG,
  ::Type{RiskParameters{T}},
  rpriors::RiskPriors{T}) where {
  T <: DiseaseStateSequence}
  sparks = Float64[rand(rng, x) for x in rpriors.sparks]
  susceptibility = Float64[rand(rng, x) for x in rpriors.susceptibility]
  infectivity = Float64[rand(rng, x) for x in rpriors.infectivity]
  transmissibility = Float64[rand(rng, x) for x in rpriors.transmissibility]
  if State_E ∈ T
    latency = Float64[rand(rng, x) for x in rpriors.latency]
  end
  if State_R ∈ T
    removal = Float64[rand(rng, x) for x in rpriors.removal]
  end
  if T == SEIR
    return RiskParameters{T}(sparks,
                             susceptibility,
                             infectivity,
                             transmissibility,
                             latency,
                             removal)
  elseif T == SEI
    return RiskParameters{T}(sparks,
                             susceptibility,
                             infectivity,
                             transmissibility,
                             latency)
  elseif T == SIR
    return RiskParameters{T}(sparks,
                             susceptibility,
                             infectivity,
                             transmissibility,
                             removal)
  elseif T == SI
    return RiskParameters{T}(sparks,
                             susceptibility,
                             infectivity,
                             transmissibility)
  end
end

function generate(rng::AbstractRNG, ::Type{RiskParameters}, mcmc::MCMC{T}) where T <: DiseaseStateSequence
  return generate(rng, RiskParameters{T}, mcmc.risk_priors)
end

function generate(
  rng::AbstractRNG,
  ::Type{RiskParameters{T}},
  last_rparams::RiskParameters{T},
  Σ::Array{Float64, 2}) where {
  T <: DiseaseStateSequence}
  rparams_vector = rand(rng, MvNormal(convert(Vector{Float64}, last_rparams), Σ))
  return _like(last_rparams, rparams_vector)
end

function generate(
  rng::AbstractRNG,
  ::Type{RiskParameters{T}},
  mc::MarkovChain{T},
  Σ::Array{Float64, 2}) where {
  T <: DiseaseStateSequence}
  return generate(rng, RiskParameters{T}, mc.risk_parameters[end], Σ)
end

function generate(
  rng::AbstractRNG,
  ::Type{SM},
  priors::Vector{UnivariateDistribution}) where {
  SM <: NucleicAcidSubstitutionModel}
  return SM([rand(rng, x) for x in priors])
end

function generate(
  rng::AbstractRNG,
  ::Type{SM},
  last_sm::SM,
  Σ::Array{Float64, 2}) where {
  SM <: NucleicAcidSubstitutionModel}
  if size(Σ) == (0,0)
    return last_sm
  else
    return SM(rand(rng, MvNormal([getproperty(last_sm, θ) for θ in [fieldnames(SM)...]], Σ)), safe=false)
  end
end

# Calls without an RNG retain the default-stream API.
generate(::Type{T}, args...; kwargs...) where {T <: Union{Transmission, TransmissionNetwork, Events, Event, RiskParameters, NucleicAcidSubstitutionModel}} =
  generate(_default_rng(), T, args...; kwargs...)
