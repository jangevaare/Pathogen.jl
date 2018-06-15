function generate(::Type{Event{T}},
                  rates::EventRates{T},
                  time::Float64=0.0) where T <: EpidemicModel
  totals = [sum(rates[state]) for state in _state_progressions[T][2:end]]
  ctotals = cumsum(totals)
  if ctotals[end] == Inf
    time = time
    new_state = _state_progressions[T][findfirst(ctotals .== Inf) + 1]
    id = findfirst(rates[new_state] .== Inf)
    return Event{T}(time, id, new_state)
  elseif ctotals[end] == 0.
    time = Inf
    return Event{T}(time)
  else
    # Generate event time
    time += rand(Exponential(1.0 / ctotals[end]))
    # Generate new state
    new_state_index = findfirst(rand(Multinomial(1, totals ./ ctotals[end])))
    new_state = _state_progressions[T][new_state_index+1]
    # Generate event individual
    id = findfirst(rand(Multinomial(1, rates[new_state] ./ totals[new_state_index])))
    return Event{T}(time, id, new_state)
  end
end

function generate(::Type{Transmission}, tr::TransmissionRates, event::Event{T}) where T <: EpidemicModel
  id = event.individual
  if _new_transmission(event)
    external_or_internal = [tr.external[id]; sum(tr.internal[:,id])]
    if rand() <= external_or_internal[1]/sum(external_or_internal)
      return ExogenousTransmission(id)
    else
      source = findfirst(rand(Multinomial(1, tr.internal[:,id] ./ external_or_internal[2])))
      return EndogenousTransmission(id, source)
    end
  else
    return NoTransmission()
  end
end

function generate(::Type{Events{T}},
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
  return generate(Events{T}, mcmc.event_observations, mcmc.event_extents)
end

function generate(::Type{Event{T}},
                  last_event::Event{T},
                  σ::Float64,
                  extents::EventExtents{T},
                  obs::EventObservations,
                  events::Events{T}) where T <: EpidemicModel
  individual = last_event.individual
  new_state = individual.new_state
  if new_state == State_E
    lb = events.infection[individual] - extents.exposure
    ub = events.infection[individual]
  elseif new_state == State_I
    if T in [SEIR; SEI]
      lb = maximum([obs.infection[individual] - extents.infection
                    events.exposure[individual]])
      ub = minimum([obs.infection[individual]
                    events.exposure[individual] + extents.exposure])
    elseif T in [SIR; SI]
      lb = obs.infection[individual] - extents.infection
      ub = obs.infection[individual]
    end
  elseif new_state == State_R
    lb = maximum([observation.removal[individual] - extents.removal
                  events.infection[individual]])
    ub = observation.removal[individual]
  end
  time = rand(TruncatedNormal(last_event.time, σ, lb, ub))
  return Event{T}(time, individual, new_state)
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
