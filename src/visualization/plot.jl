@recipe function plot(events::Events{T}, min::Float64, max::Float64) where T<: EpidemicModel
  xguide --> "Time"
  yguide --> "N"
  xlims --> (min - 1.0, max + 1.0)
  @series begin
    seriestype := :steppost
    seriescolor --> :purple
    time, count = _epidemic_curve(events, State_S, min, max)
    label := "S"
    time, count
  end
  if T ∈ [SEIR; SEI]
    @series begin
      seriestype := :steppost
      seriescolor --> :lightblue4
      time, count = _epidemic_curve(events, State_E, min, max)
      label := "E"
      time, count
    end
  end
  @series begin
    seriestype := :steppost
    seriescolor --> :lightgreen
    time, count = _epidemic_curve(events, State_I, min, max)
    label := "I"
    time, count
  end
  if T ∈ [SEIR; SIR]
    @series begin
      seriestype := :steppost
      seriescolor --> :yellow
      time, count = _epidemic_curve(events, State_R, min, max)
      label := "R"
      time, count
    end
  end
end

function plot(events::Events{T}) where T <: EpidemicModel
  return plot(events, minimum(events), maximum(events))
end

@recipe function plot(network::TransmissionNetwork, pop::DataFrame, events::Events{T},  time::Float64) where T <: EpidemicModel
  xguide --> ""
  yguide --> ""
  legend --> :topright
  xlims --> extrema(pop[:x]) .+ (sum(extrema(pop[:x]).*(-1,1)) .* (-0.05, 0.05))
  ylims --> extrema(pop[:y]) .+ (sum(extrema(pop[:y]).*(-1,1)) .* (-0.05, 0.05))
  ids_susceptible = _ids_by_state(events, State_S, time)
  if T ∈ [SEIR; SEI]
    ids_exposed = _ids_by_state(events, State_E, time)
    for i ∈ ids_exposed
      @series begin
        x, y = _pathway_to(pop, network, i)
        seriestype := :path
        seriescolor := :black
        label := ""
        x, y
      end
    end
  end
  ids_infected = _ids_by_state(events, State_I, time)
  for i ∈ ids_infected
    @series begin
      x, y = _pathway_to(pop, network, i)
      seriestype := :path
      seriescolor := :black
      label := ""
      x, y
    end
  end
  if T ∈ [SEIR; SIR]
    ids_removed = _ids_by_state(events, State_R, time)
    for i ∈ ids_removed
      @series begin
        x, y = _pathway_to(pop, network, i)
        seriestype := :path
        seriescolor := :black
        label := ""
        x, y
      end
    end
  end
  @series begin
    x, y = _population_plot(pop, ids_susceptible)
    seriestype := :scatter
    seriescolor --> :purple
    label := "S"
    x, y
  end
  if T ∈ [SEIR; SEI]
    @series begin
      x, y = _population_plot(pop, ids_exposed)
      seriestype := :scatter
      seriescolor --> :lightblue4
      label := "E"
      x, y
    end
  end
  @series begin
    x, y = _population_plot(pop, ids_infected)
    seriestype := :scatter
    seriescolor --> :lightgreen
    label := "I"
    x, y
  end
  if T ∈ [SEIR; SIR]
    @series begin
      x, y = _population_plot(pop, ids_removed)
      seriestype := :scatter
      seriescolor --> :yellow
      label := "R"
      x, y
    end
  end
end
