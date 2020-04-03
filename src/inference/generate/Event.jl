function _pathway_to(id::Int64,
                     network::TransmissionNetwork;
                     depth::Real=Inf)
  path = Int64[]
  if network.external[id] | any(network.internal[:, id])
    push!(path, id)
    next_id = findfirst(network.internal[:, id])
    while (length(path) <= depth)&& (typeof(next_id) != Nothing)
      push!(path, next_id)
      next_id = findfirst(network.internal[:, path[end]])
    end
  end
  @debug "_pathway_to: transmission pathway to $id: $path"
  return path
end

function _bounds(id::Int64,
                 new_state::DiseaseState,
                 extents::EventExtents{S},
                 obs::EventObservations{S, M},
                 events::Events{S},
                 network::TransmissionNetwork) where {
                   S <: DiseaseStateSequence,
                   M <: ILM}
  if new_state == State_E
    lowerbounds = [events.infection[id] - extents.exposure]
    upperbounds = [events.infection[id]]
    path_to = _pathway_to(id, network, depth = 1)
    if length(path_to) > 1
      parent_host = path_to[2]
      push!(lowerbounds, events.infection[parent_host])
      if (S == SEIR)&& !isnan(events.removal[parent_host])
        push!(upperbounds, events.removal[parent_host])
      end
    end
  elseif new_state == State_I
    path_from = _pathway_from(id, network, depth = 1)
    if S in [SEIR; SEI]
      lowerbounds = [obs.infection[id] - extents.infection
                     events.exposure[id]]
      upperbounds = [obs.infection[id]
                     events.exposure[id] + extents.exposure]
      if length(path_from) > 1
        child_hosts = path_from[2:end]
        append!(upperbounds, events.exposure[child_hosts])
      end
    elseif S in [SIR; SI]
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
        if (S == SIR)&& !isnan(events.removal[parent_host])
          push!(upperbounds, events.removal[parent_host])
        end
      end
    end
  elseif new_state == State_R
    path_from = _pathway_from(id, network, depth = 1)
    lowerbounds = [obs.removal[id] - extents.removal
                   obs.infection[id]]
    upperbounds = [obs.removal[id]]
    if length(path_from) > 1
      child_hosts = path_from[2:end]
      if S == SEIR
        append!(lowerbounds, events.exposure[child_hosts])
      elseif S == SIR
        append!(lowerbounds, events.infection[child_hosts])
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
    lowerbounds = [events.infection[id] - extents.exposure]
    upperbounds = [events.infection[id]]
  elseif new_state == State_I
    if S in [SEIR; SEI]
      lowerbounds = [obs.infection[id] - extents.infection
                     events.exposure[id]]
      upperbounds = [obs.infection[id]
                     events.exposure[id] + extents.exposure]
    elseif S in [SIR; SI]
      lowerbounds = [obs.infection[id] - extents.infection]
      upperbounds = [obs.infection[id]]
    end
  elseif new_state == State_R
    lowerbounds = [obs.removal[id] - extents.removal
                   obs.infection[id]]
    upperbounds = [obs.removal[id]]
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
