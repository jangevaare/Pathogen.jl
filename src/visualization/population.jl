function _population(pop::Population, ids)
  x = Float64[]
  y = Float64[]
  for i in ids
    push!(x, pop.risks[i, :x])
    push!(y, pop.risks[i, :y])
  end
  return x, y
end

function _ids_by_state(events::Events{T},
                       state::DiseaseState,
                       time::Float64) where T <: DiseaseStateSequence
  if time < 0.0
    @error "Time must be ≥ 0.0"
  elseif state ∉ convert(DiseaseStates, T)
    @error "Invalid state specified"
  end
  ids = Int64[]
  if state == convert(DiseaseStates, T)[1] # S
    nextstate = advance(state, T) # Either E or I
    append!(ids, findall(events[nextstate] .> Ref(time))) # E/I after `time`
    append!(ids, findall(isnan.(events[nextstate]))) # Never E/I
  elseif state in convert(DiseaseStates, T)[2:end-1]
    nextstate = advance(state, T) # Either I or R
    append!(ids, findall((events[state] .<= Ref(time)) .& (events[nextstate] .> Ref(time)))) # E/I at or before time and I/R after time
    append!(ids, findall((events[state] .<= Ref(time)) .& isnan.(events[nextstate]))) # E/I at or before time and never I/R
  elseif state == convert(DiseaseStates, T)[end] # I or R
    append!(ids, findall(events[state] .<= Ref(time))) # I/R at or before time
  end
  @debug "Individual(s) in state $state at t = $time" Individuals = ids
  return ids
end

@recipe function f(pop::Population) where T <: DiseaseStateSequence
  xguide --> ""
  yguide --> ""
  legend --> :none
  xlims --> extrema(pop.risks[!, :x]) .+ (sum(extrema(pop.risks[!, :x]).*(-1,1)) .* (-0.05, 0.05))
  ylims --> extrema(pop.risks[!, :y]) .+ (sum(extrema(pop.risks[!, :y]).*(-1,1)) .* (-0.05, 0.05))
  aspect_ratio := :equal
  seriestype := :scatter
  markerstrokecolor --> :black
  markercolor --> :black
  markersize --> 2.75
  axis --> nothing
  titlefontcolor --> :black
  framestyle --> :none
  _population(pop, 1:individuals(pop))
end


@recipe function f(pop::Population,
                   events::Events{T},
                   time::Float64) where T <: DiseaseStateSequence
  xguide --> ""
  yguide --> ""
  legend --> :topright
  xlims --> extrema(pop.risks[!, :x]) .+ (sum(extrema(pop.risks[!, :x]).*(-1,1)) .* (-0.05, 0.05))
  ylims --> extrema(pop.risks[!, :y]) .+ (sum(extrema(pop.risks[!, :y]).*(-1,1)) .* (-0.05, 0.05))
  aspect_ratio := :equal
  markerstrokecolor --> :black
  markersize --> 2.75
  axis --> nothing
  titlefontcolor --> :black
  framestyle --> :none
  @series begin
    ids_susceptible = _ids_by_state(events, State_S, time)
    x, y = _population(pop, ids_susceptible)
    seriestype := :scatter
    seriescolor --> :purple
    label --> "S"
    x, y
  end
  if T in [SEIR; SEI]
    @series begin
      ids_exposed = _ids_by_state(events, State_E, time)
      x, y = _population(pop, ids_exposed)
      seriestype := :scatter
      seriescolor --> :lightblue4
      label --> "E"
      x, y
    end
  end
  @series begin
    ids_infected = _ids_by_state(events, State_I, time)
    x, y = _population(pop, ids_infected)
    seriestype := :scatter
    seriescolor --> :lightgreen
    label --> "I"
    x, y
  end
  if T in [SEIR; SIR]
    @series begin
      ids_removed = _ids_by_state(events, State_R, time)
      x, y = _population(pop, ids_removed)
      seriestype := :scatter
      seriescolor --> :yellow
      label --> "R"
      x, y
    end
  end
end
