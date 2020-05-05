function advance(x::DiseaseState, ::Type{S}) where {S <: DiseaseStateSequence}
  states = convert(DiseaseStates, S)
  current_index = findfirst(Ref(x) .== states)
  return states[current_index + 1]
end

function _count_by_state(events::Events{T},
                         state::DiseaseState,
                         time::Float64) where T <: DiseaseStateSequence
  local n_ids
  if state ∉ convert(DiseaseStates, T)
    @error "Invalid state specified"
  end
  if state == convert(DiseaseStates, T)[1] # S
    nextstate = advance(state, T) # Either E or I
    n_ids = sum(events[nextstate] .> Ref(time)) # E/I after `time`
    n_ids += sum(isnan.(events[nextstate])) # Never E/I
  elseif state in convert(DiseaseStates, T)[2:end-1]
    nextstate = advance(state, T) # Either I or R
    n_ids = sum((events[state] .<= Ref(time)) .& (events[nextstate] .> Ref(time))) # E/I at or before time and I/R after time
    n_ids += sum((events[state] .<= Ref(time)) .& isnan.(events[nextstate])) # E/I at or before time and never I/R
  elseif state == convert(DiseaseStates, T)[end] # I or R
    n_ids = sum(events[state] .<= Ref(time)) # I/R at or before time
  end
  @debug "$n_ids individual(s) in state $state at t = $time"
  return n_ids
end

function _epidemic_curve(events::Events{T},
                         state::DiseaseState,
                         min::Float64,
                         max::Float64) where T <: DiseaseStateSequence
  if min >= max
    @error "Minimum time must be less than maximum time"
  end
  local times
  if state == convert(DiseaseStates, T)[1]
    nextstate = convert(DiseaseStates, T)[2]
    times = events[nextstate]
  elseif state in convert(DiseaseStates, T)[2:end-1]
    nextstate = advance(state, T)
    times = events[[state; nextstate]][:]
  elseif state == convert(DiseaseStates, T)[end]
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

@recipe function f(events::Events{T},
                   state::DiseaseState,
                   min::Float64,
                   max::Float64) where T<: DiseaseStateSequence
  xguide --> "Time"
  yguide --> "N"
  xlims --> (min - 1.0, max + 1.0)
  linewidth --> 2.0
  linecolor --> :cornflowerblue
  label --> ""
  seriestype --> :steppost
  _epidemic_curve(events, state, min, max)
end

@recipe function f(events::Events{T},
                   state::DiseaseState) where T <: DiseaseStateSequence
  events, state, 0.0, maximum(events)
end

@recipe function f(events::Events{T},
                   min::Float64,
                   max::Float64) where T<: DiseaseStateSequence
  @series begin
    linecolor --> :purple
    label --> "S"
    events, State_S, min, max
  end
  if T in [SEIR; SEI]
    @series begin
      linecolor --> :lightblue4
      label --> "E"
      events, State_E, min, max
    end
  end
  @series begin
    linecolor --> :lightgreen
    label --> "I"
    events, State_I, min, max
  end
  if T in [SEIR; SIR]
    @series begin
      linecolor --> :yellow
      label --> "R"
      events, State_R, min, max
    end
  end
end

@recipe function f(events::Events{T}) where T <: DiseaseStateSequence
  events, 0.0, maximum(events)
end
