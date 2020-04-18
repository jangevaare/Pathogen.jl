function _bounds(id::I,
                 new_state::DiseaseState,
                 extents::EventExtents{T},
                 obs::EventObservations{T},
                 events::Events{T},
                 network::TransmissionNetwork) where {I <: Integer, T <: EpidemicModel}
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
                   obs.infection[id]]
    upperbounds = [obs.removal[id]]
    if length(path_from) > 1
      child_hosts = path_from[2:end]
      if T == SEIR
        append!(lowerbounds, events.exposure[child_hosts])
      elseif T == SIR
        append!(lowerbounds, events.infection[child_hosts])
      end
    end
  end
  @debug "Transition of $id into $new_state with bounds: \n  [max($(round.(lowerbounds, digits=3))),\n   min($(round.(upperbounds, digits=3)))]"
  lowerbound = maximum(lowerbounds)
  upperbound = minimum(upperbounds)
  if lowerbound >= upperbound
    @warn "Lower >= upper bound for transition of i=$id into $new_state" lower = lowerbounds upper = upperbounds
  end
  return lowerbound, upperbound
end

function _bounds(last_event::Event{T},
                 extents::EventExtents{T},
                 obs::EventObservations{T},
                 events::Events{T},
                 network::TransmissionNetwork) where T <: EpidemicModel
  return _bounds(last_event.individual,
                 last_event.new_state,
                 extents, obs, events, network)
end


function _bounds(id::I,
                 new_state::DiseaseState,
                 extents::EventExtents{T},
                 obs::EventObservations{T},
                 events::Events{T}) where {I <: Integer, T <: EpidemicModel}
  if new_state == State_E
    lowerbounds = [events.infection[id] - extents.exposure]
    upperbounds = [events.infection[id]]
  elseif new_state == State_I
    if T in [SEIR; SEI]
      lowerbounds = [obs.infection[id] - extents.infection
                     events.exposure[id]]
      upperbounds = [obs.infection[id]
                     events.exposure[id] + extents.exposure]
    elseif T in [SIR; SI]
      lowerbounds = [obs.infection[id] - extents.infection]
      upperbounds = [obs.infection[id]]
    end
  elseif new_state == State_R
    lowerbounds = [obs.removal[id] - extents.removal
                   obs.infection[id]]
    upperbounds = [obs.removal[id]]
  end
  @debug "Transition of $id into $new_state with bounds: \n  [max($(round.(lowerbounds, digits=3))),\n   min($(round.(upperbounds, digits=3)))]"
  lowerbound = maximum(lowerbounds)
  upperbound = minimum(upperbounds)
  if lowerbound >= upperbound
    @warn "Lower >= upper bound for transition of i=$id into $new_state" lower = lowerbounds upper = upperbounds
  end
  return lowerbound, upperbound
end

function _bounds(last_event::Event{T},
                 extents::EventExtents{T},
                 obs::EventObservations{T},
                 events::Events{T}) where T <: EpidemicModel
  return _bounds(last_event.individual,
                 last_event.new_state,
                 extents, obs, events)
end
