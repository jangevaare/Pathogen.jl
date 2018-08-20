function _epidemic_curve(events::Events{T}, state::DiseaseState, min::Float64, max::Float64) where T <: EpidemicModel
  if min >= max
    @error "Minimum time must be less than maximum time"
  end
  local times
  if state == _state_progressions[T][1]
    nextstate = _state_progressions[T][2]
    times = events[nextstate]
  elseif state in _state_progressions[T][2:end-1]
    nextstate = advance(state, T)
    times = events[[state; nextstate]][:]
  elseif state == _state_progressions[T][end]
    times = events[state]
  else
    @error "Invalid state specified"
  end
  times = times[Ref(min) .< times .< Ref(max)]
  sort!(times)
  insert!(times, 1, min)
  push!(times, max)
  counts = [_count_by_state(events, state, t) for t in times]
  return times, counts
end

function _epidemic_curve(events::Events{T}, state::DiseaseState) where T <: EpidemicModel
  return _epidemic_curve(events, state, 0.0, maximum(events))
end
