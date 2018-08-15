function generate(::Type{Event},
                  rates::EventRates{T},
                  time::Float64=0.0) where T <: EpidemicModel
  totals = Weights([sum(rates[state]) for state in _state_progressions[T][2:end]])
  if sum(totals) == Inf
    new_state = sample(_state_progressions[T][2:end], totals)
    id = sample(1:rates.individuals, Weights(rates[new_state]))
    return Event{T}(time, id, new_state)
  elseif sum(totals) == 0.0
    return Event{T}(Inf)
  else
    # Generate new state
    new_state = sample(_state_progressions[T][2:end], totals)
    # Generate event individual
    id = sample(1:rates.individuals, Weights(rates[new_state]))
    return Event{T}(time + rand(Exponential(1.0 / sum(totals))), id, new_state)
  end
end

function generate(::Type{Transmission},
                  tr::TransmissionRates,
                  event::Event{T}) where T <: EpidemicModel
  id = event.individual
  if _new_transmission(event)
    external_or_internal = Weights([tr.external[id]; sum(tr.internal[:,id])])
    if sum(external_or_internal) == 0.0
      @debug "All transmission rates = 0.0, exogenous transmission generated"
      return ExogenousTransmission(id)
    elseif sample([true; false], external_or_internal)
      @debug "Exogenous tranmission generated"
      return ExogenousTransmission(id)
    else
      source = sample(1:tr.individuals, Weights(tr.internal[:, id]))
      @debug "Endogenous transmission generated (source id = $source)"
      return EndogenousTransmission(id, source)
    end
  else
    @logmsg LogLevel(-3000) "No transmission generated"
    return NoTransmission()
  end
end

function generate(::Type{Events},
                  obs::EventObservations{T},
                  extents::EventExtents{T}) where T <: EpidemicModel
  events = Events{T}(obs.individuals)
  exposed_state = T in [SEIR; SEI]
  removed_state = T in [SEIR; SIR]
  for i = 1:obs.individuals
    if !isnan(obs.infection[i])
      i_lb = obs.infection[i] - extents.infection
      i_ub = obs.infection[i]
      i_time = rand(Uniform(i_lb, i_ub))
      update!(events, Event{T}(i_time, i, State_I))
      if exposed_state
        e_lb = i_time - extents.exposure
        e_ub = i_time
        e_time = rand(Uniform(e_lb, e_ub))
        update!(events, Event{T}(e_time, i, State_E))
      end
      if removed_state && !isnan(obs.removal[i])
        r_lb = maximum([obs.infection[i]; obs.removal[i] - extents.removal])
        r_ub = obs.removal[i]
        r_time = rand(Uniform(r_lb, r_ub))
        update!(events, Event{T}(r_time, i, State_R))
      end
    end
  end
  return events
end

function generate(::Type{Events}, mcmc::MCMC{T}) where T <: EpidemicModel
  return generate(Events, mcmc.event_observations, mcmc.event_extents)
end

function generate(::Type{Event},
                  last_event::Event{T},
                  σ::Float64,
                  extents::EventExtents{T},
                  obs::EventObservations{T},
                  events::Events{T}) where T <: EpidemicModel
  id = last_event.individual
  new_state = last_event.new_state
  if new_state == State_E
    lb = events.infection[id] - extents.exposure
    ub = events.infection[id]
  elseif new_state == State_I
    if T in [SEIR; SEI]
      lb = maximum([obs.infection[id] - extents.infection
                    events.exposure[id]])
      ub = minimum([obs.infection[id]
                    events.exposure[id] + extents.exposure])
    elseif T in [SIR; SI]
      lb = obs.infection[id] - extents.infection
      ub = obs.infection[id]
    end
  elseif new_state == State_R
    lb = maximum([obs.removal[id] - extents.removal
                  events.infection[id]])
    ub = obs.removal[id]
  end
  time = rand(TruncatedNormal(last_event.time, σ, lb, ub))
  return Event{T}(time, id, new_state)
end

function generate(::Type{Event},
                  last_event::Event{T},
                  σ::Float64,
                  extents::EventExtents{T},
                  obs::EventObservations,
                  events::Events{T},
                  network::TransmissionNetwork) where T <: EpidemicModel
  lowerbound, upperbound = _bounds(last_event, extents, obs, events, network)
  time = rand(TruncatedNormal(last_event.time,
                              σ,
                              lowerbound,
                              upperbound))
  return Event(time, last_event)
end

function generate(::Type{RiskParameters{T}}, rpriors::RiskPriors{T}) where T <: EpidemicModel
  sparks = Float64[rand(x) for x in rpriors.sparks]
  susceptibility = Float64[rand(x) for x in rpriors.susceptibility]
  transmissibility = Float64[rand(x) for x in rpriors.transmissibility]
  infectivity = Float64[rand(x) for x in rpriors.infectivity]
  if T in [SEIR; SEI]
    latency = Float64[rand(x) for x in rpriors.latency]
  end
  if T in [SEIR; SIR]
    removal = Float64[rand(x) for x in rpriors.removal]
  end
  if T == SEIR
    return RiskParameters{T}(sparks,
                             susceptibility,
                             transmissibility,
                             infectivity,
                             latency,
                             removal)
  elseif T == SEI
    return RiskParameters{T}(sparks,
                             susceptibility,
                             transmissibility,
                             infectivity,
                             latency)
  elseif T == SIR
    return RiskParameters{T}(sparks,
                             susceptibility,
                             transmissibility,
                             infectivity,
                             removal)
  elseif T == SI
    return RiskParameters{T}(sparks,
                             susceptibility,
                             transmissibility,
                             infectivity)
  end
end

function generate(::Type{RiskParameters}, mcmc::MCMC{T}) where T <: EpidemicModel
  return generate(RiskParameters{T}, mcmc.risk_priors)
end

function generate(::Type{RiskParameters{T}},
                  last_rparams::RiskParameters{T},
                  Σ::Array{Float64, 2}) where T <: EpidemicModel
  rparams_vector = rand(MvNormal(convert(Vector, last_rparams), Σ))
  return _like(last_rparams, rparams_vector)
end

function generate(::Type{RiskParameters{T}},
                  mc::MarkovChain{T},
                  Σ::Array{Float64, 2}) where T <: EpidemicModel
  return generate(RiskParameters{T}, mc.risk_parameters[end], Σ)
end
