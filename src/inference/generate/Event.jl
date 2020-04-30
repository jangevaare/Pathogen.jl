function _bounds(id::Int64,
                 new_state::DiseaseState,
                 extents::EventExtents{S},
                 obs::EventObservations{S, M},
                 events::Events{S},
                 network::TransmissionNetwork) where {
                 S <: DiseaseStateSequence,
                 M <: ILM}
  if new_state == State_E
    lowerbounds = [events.infection[id] - extents.exposure[2]]
    upperbounds = [events.infection[id] - extents.exposure[1]]
    source = _source(id, network)
    if source !== nothing
      push!(lowerbounds, events.infection[source])
      if (S == SEIR) && !isnan(events.removal[source])
        push!(upperbounds, events.removal[source])
      end
    end
  elseif new_state == State_I
    path_from = _pathway_from(id, network, depth = 1, include_id = false)
    if S in [SEIR; SEI]
      lowerbounds = [obs.infection[id] - extents.infection[2]
                     events.exposure[id] + extents.exposure[1]]
      upperbounds = [obs.infection[id] - extents.infection[1]
                     events.exposure[id] + extents.exposure[2]]
      if length(path_from) > 0
        append!(upperbounds, events.exposure[path_from])
      end
    elseif S in [SIR; SI]
      lowerbounds = [obs.infection[id] - extents.infection[2]]
      upperbounds = [obs.infection[id] - extents.infection[1]]
      if length(path_from) > 0
        append!(upperbounds, events.infection[path_from])
      end
      source = _source(id, network)
      if source !== nothing
        push!(lowerbounds, events.infection[source])
        if (S == SIR) && !isnan(events.removal[source])
          push!(upperbounds, events.removal[source])
        end
      end
    end
  elseif new_state == State_R
    path_from = _pathway_from(id, network, depth = 1, include_id = false)
    lowerbounds = [obs.removal[id] - extents.removal[2]
                   obs.infection[id]]
    upperbounds = [obs.removal[id] - extents.removal[1]]
    if length(path_from) > 0
      if S == SEIR
        append!(lowerbounds, events.exposure[path_from])
      elseif S == SIR
        append!(lowerbounds, events.infection[path_from])
      end
    end
  end
  # Must be later than epidemic start time
  push!(lowerbounds, obs.start_time)
  @debug "Transition of $id into $new_state with bounds: \n  [max($(round.(lowerbounds, digits=3))),\n   min($(round.(upperbounds, digits=3)))]"
  lowerbound = maximum(lowerbounds)
  upperbound = minimum(upperbounds)
  if lowerbound >= upperbound
    @warn "Lower >= upper bound for transition of i=$id into $new_state" lower = lowerbounds upper = upperbounds
  end
  return lowerbound, upperbound
end

function _bounds(last_event::Event{S},
                 extents::EventExtents{S},
                 obs::EventObservations{S, M},
                 events::Events{S},
                 network::TransmissionNetwork) where {
                 S <: DiseaseStateSequence,
                 M <: ILM}
  return _bounds(last_event.individual,
                 last_event.new_state,
                 extents, obs, events, network)
end


function _bounds(id::Int64,
                 new_state::DiseaseState,
                 extents::EventExtents{S},
                 obs::EventObservations{S, M},
                 events::Events{S}) where {
                 S <: DiseaseStateSequence,
                 M <: ILM}
  if new_state == State_E
    lowerbounds = [events.infection[id] - extents.exposure[2]]
    upperbounds = [events.infection[id] - extents.exposure[1]]
  elseif new_state == State_I
    if S in [SEIR; SEI]
      lowerbounds = [obs.infection[id] - extents.infection[2]
                     events.exposure[id] + extents.exposure[1]]
      upperbounds = [obs.infection[id] - extents.infection[1]
                     events.exposure[id] + extents.exposure[2]]
    elseif S in [SIR; SI]
      lowerbounds = [obs.infection[id] - extents.infection[2]]
      upperbounds = [obs.infection[id] - extents.infection[1]]
    end
  elseif new_state == State_R
    lowerbounds = [obs.removal[id] - extents.removal[2]
                   obs.infection[id]]
    upperbounds = [obs.removal[id] - extents.removal[1]]
  end
  # Must be later than epidemic start time
  push!(lowerbounds, obs.start_time)
  @debug "Transition of $id into $new_state with bounds: \n  [max($(round.(lowerbounds, digits=3))),\n   min($(round.(upperbounds, digits=3)))]"
  lowerbound = maximum(lowerbounds)
  upperbound = minimum(upperbounds)
  if lowerbound >= upperbound
    @warn "Lower >= upper bound for transition of i=$id into $new_state" lower = lowerbounds upper = upperbounds
  end
  return lowerbound, upperbound
end

function _bounds(last_event::Event{S},
                 extents::EventExtents{S},
                 obs::EventObservations{S, M},
                 events::Events{S}) where {
                 S <: DiseaseStateSequence,
                 M <: ILM}
  return _bounds(last_event.individual,
                 last_event.new_state,
                 extents, obs, events)
end

function generate(::Type{Event},
                  last_event::Event{S},
                  σ::Float64,
                  extents::EventExtents{S},
                  obs::EventObservations{S, M},
                  events::Events{S}) where {
                  S <: DiseaseStateSequence,
                  M <: ILM}
  lowerbound, upperbound = _bounds(last_event, extents, obs, events)
  time = rand(truncated(Normal(_time(last_event), σ),
                        lowerbound,
                        upperbound))
  return Event(time, last_event)
end

function generate(::Type{Event},
                  last_event::Event{S},
                  σ::Float64,
                  extents::EventExtents{S},
                  obs::EventObservations,
                  events::Events{S},
                  network::TransmissionNetwork) where {
                  S <: DiseaseStateSequence}
  lowerbound, upperbound = _bounds(last_event, extents, obs, events, network)
  time = rand(truncated(Normal(_time(last_event), σ),
                        lowerbound,
                        upperbound))
  return Event(time, last_event)
end
