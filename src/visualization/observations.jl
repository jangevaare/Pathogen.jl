function _count_by_state(events::EventObservations{S, M},
                         state::DiseaseState,
                         time::Float64) where {
                         S <: DiseaseStateSequence,
                         M <: ILM}
  local n_ids
  if state ∉ [State_I, State_R]
    throw(ErrorException("Invalid state specified"))
  elseif state == State_I && S ∈ [SEIR, SIR]
    nextstate = State_R
    n_ids = sum((events[state] .<= Ref(time)) .& (events[nextstate] .> Ref(time))) # E/I at or before time and I/R after time
    n_ids += sum((events[state] .<= Ref(time)) .& isnan.(events[nextstate])) # E/I at or before time and never I/R 
  else
    n_ids = sum(events[state] .<= Ref(time)) # I/R at or before time
  end
  @debug "$n_ids individual(s) observered in state $state at t = $time"
  return n_ids
end

function _obs_curve(events::EventObservations{S, M},
                    state::DiseaseState,
                    min::Float64,
                    max::Float64) where {
                    S <: DiseaseStateSequence,
                    M <: ILM}
  if min >= max
    error("Minimum time must be less than maximum time")
  end
  local times
  if state == State_I && T ∈ [SEIR, SIR]
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

@recipe function f(events::EventObservations{S, M},
                   state::DiseaseState,
                   min::Float64,
                   max::Float64) where {
                   S <: DiseaseStateSequence,
                   M <: ILM}
  xguide --> "Time"
  yguide --> "N"
  xlims --> (min - 1.0, max + 1.0)
  linewidth --> 2.0
  linecolor --> :cornflowerblue
  label --> ""
  seriestype --> :steppost
  _obs_curve(events, state, min, max)
 end

@recipe function f(events::EventObservations{S, M},
                   state::DiseaseState) where {
                   S <: DiseaseStateSequence,
                   M <: ILM}
   events, state, 0.0, maximum(events)
end

@recipe function f(events::EventObservations{S, M},
                    min::Float64,
                    max::Float64) where {
                    S <: DiseaseStateSequence,
                    M <: ILM}
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

@recipe function f(events::EventObservations{S, M}) where {
                   S <: DiseaseStateSequence,
                   M <: ILM}
  events, 0.0, maximum(events)
end