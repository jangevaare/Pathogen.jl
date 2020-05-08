function _count_by_state(events::EventObservations{T},
                         state::DiseaseState,
                         time::Float64) where T <: DiseaseStateSequence
  local n_ids
  if state == State_I && State_R ∈ T
    nextstate = State_R
    n_ids = sum((events[state] .<= Ref(time)) .& (events[nextstate] .> Ref(time))) # E/I at or before time and I/R after time
    n_ids += sum((events[state] .<= Ref(time)) .& isnan.(events[nextstate])) # E/I at or before time and never I/R 
  else
    n_ids = sum(events[state] .<= Ref(time)) # I/R at or before time
  end
  @debug "$n_ids individual(s) observered in state $state at t = $time"
  return n_ids
end

function _obs_curve(events::EventObservations{T},
                         state::DiseaseState,
                         min::Float64,
                         max::Float64) where T <: DiseaseStateSequence
  if min >= max
    @error "Minimum time must be less than maximum time"
  end
  local times
  if state == State_I && State_R ∈ T
    times = events[[State_I; State_R]][:]
  else
    times = events[state]
  end
  times = times[Ref(min) .< times .< Ref(max)]
  sort!(times)
  insert!(times, 1, min)
  push!(times, max)
  counts = [_count_by_state(events, state, t) for t in times]
  return times, counts
end

@recipe function f(events::EventObservations{T},
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
  _obs_curve(events, state, min, max)
end

@recipe function f(events::EventObservations{T},
                   state::DiseaseState) where T <: DiseaseStateSequence
  events, state, 0.0, maximum(events)
end

@recipe function f(events::EventObservations{T},
                   min::Float64,
                   max::Float64) where T<: DiseaseStateSequence
  @series begin
    linecolor --> :lightgreen
    label --> "I"
    events, State_I, min, max
  end
  if State_R ∈ T
    @series begin
      linecolor --> :yellow
      label --> "R"
      events, State_R, min, max
    end
  end
end

@recipe function f(events::EventObservations{T}) where T <: DiseaseStateSequence
  events, 0.0, maximum(events)
end
