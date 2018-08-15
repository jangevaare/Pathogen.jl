function _ids_by_state(events::Events{T}, state::DiseaseState, time::Float64) where T <: EpidemicModel
  if time < 0.0
    @error "Time must be ≥ 0.0"
  elseif state ∉ _state_progressions[T]
    @error "Invalid state specified"
  end
  ids = Int64[]
  if state == _state_progressions[T][1] # S
    nextstate = advance(state, T) # Either E or I
    append!(ids, findall(events[nextstate] .> Ref(time))) # E/I after `time`
    append!(ids, findall(isnan.(events[nextstate]))) # Never E/I
  elseif state ∈ _state_progressions[T][2:end-1]
    nextstate = advance(state, T) # Either I or R
    append!(ids, findall((events[state] .<= Ref(time)) .& (events[nextstate] .> Ref(time)))) # E/I at or before time and I/R after time
    append!(ids, findall((events[state] .<= Ref(time)) .& isnan.(events[nextstate]))) # E/I at or before time and never I/R
  elseif state == _state_progressions[T][end] # I or R
    append!(ids, findall(events[state] .<= Ref(time))) # I/R at or before time
  end
  @debug "Individual(s) in state $state at t = $time" Individuals = ids
  return ids
end
