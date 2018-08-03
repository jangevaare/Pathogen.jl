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

function generate(::Type{Transmission},
                  tr::TransmissionRates,
                  event::Event{T};
                  debug_level::Int64=0) where T <: EpidemicModel
  id = event.individual
  if _new_transmission(event)
    external_or_internal = [tr.external[id]; sum(tr.internal[:,id])]
    if !any(external_or_internal .> 0.0)
      if debug_level >= 4
        println("generate: all transmission rates = 0.0, exogenous tranmission generated")
      end
      return ExogenousTransmission(id)
    elseif rand() <= external_or_internal[1]/sum(external_or_internal)
      if debug_level >= 4
        println("generate: exogenous tranmission generated")
      end
      return ExogenousTransmission(id)
    else
      source = findfirst(rand(Multinomial(1, tr.internal[:,id] ./ external_or_internal[2])))
      if debug_level >= 4
        println("generate: endogenous tranmission generated (source id = $source)")
      end
      return EndogenousTransmission(id, source)
    end
  else
    if debug_level >= 4
      println("generate: no tranmission generated")
    end
    return NoTransmission()
  end
end

function generate(::Type{Events{T}},
                  obs::EventObservations{T},
                  extents::EventExtents{T};
                  debug_level::Int64=0) where T <: EpidemicModel
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

function generate(::Type{Events}, mcmc::MCMC{T}; debug_level::Int64=0) where T <: EpidemicModel
  return generate(Events{T}, mcmc.event_observations, mcmc.event_extents, debug_level=debug_level)
end

function generate(::Type{Event{T}},
                  last_event::Event{T},
                  σ::Float64,
                  extents::EventExtents{T},
                  obs::EventObservations,
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

function generate(::Type{Event{T}},
                  last_event::Event{T},
                  σ::Float64,
                  extents::EventExtents{T},
                  obs::EventObservations,
                  events::Events{T},
                  network::TransmissionNetwork;
                  debug_level::Int64 = 0) where T <: EpidemicModel
  id = last_event.individual
  new_state = last_event.new_state
  if new_state == State_E
    lowerbounds = [events.infection[id] - extents.exposure]
    upperbounds = [events.infection[id]]
    path_to = _pathway_to(id, network, depth = 1)
    if length(path_to) > 1
      parent_host = path_to[2]
      push!(lowerbounds, events.infection[parent_host])
      if (T == SEIR) && !isnan(events.removal[parent_host])
        push!(upperbounds, events.removal[parent_host])
      end
    end
  elseif new_state == State_I
    path_from = _pathway_from(id, network, depth = 1)
    if T in [SEIR; SEI]
      lowerbounds = [obs.infection[id] - extents.infection
                     events.exposure[id]]
      upperbounds = [obs.infection[id]
                     events.exposure[id] + extents.exposure]
      if length(path_from) > 1
        child_hosts = path_from[2:end]
        append!(upperbounds, events.exposure[child_hosts])
      end
    elseif T in [SIR; SI]
      lowerbounds = [obs.infection[id] - extents.infection]
      upperbounds = [obs.infection[id]]
      path_to = _pathway_to(id, network, depth = 1)
      if length(path_from) > 1
        child_hosts = path_from[2:end]
        append!(upperbounds, events.infection[child_hosts])
      end
      if length(path_to) > 1
        parent_host = path_to[2]
        push!(lowerbounds, events.infection[parent_host])
        if (T == SIR) && !isnan(events.removal[parent_host])
          push!(upperbounds, events.removal[parent_host])
        end
      end
    end
  elseif new_state == State_R
    path_from = _pathway_from(id, network, depth = 1)
    lowerbounds = [obs.removal[id] - extents.removal
                   events.infection[id]]
    upperbounds = [obs.removal[id]]
    if length(path_from) > 1
      child_hosts = path_from[2:end]
      if T == SEIR
        append!(lowerbounds, events.exposure[child_hosts])
      elseif T == SIR
        append!(lowerbounds, events.exposure[child_hosts])
      end
    end
  end
  if debug_level >= 4
    println("Transition of $id into $new_state with bounds: \n[max($(round(lowerbounds, 3))),\n min($(round(upperbounds, 3)))]")
  end
  time = rand(TruncatedNormal(last_event.time,
                              σ,
                              maximum(lowerbounds),
                              minimum(upperbounds)))
  return Event{T}(time, id, new_state)
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
