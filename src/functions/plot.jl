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
  if T in [SEIR; SEI]
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
  if T in [SEIR; SIR]
    @series begin
      seriestype := :steppost
      seriescolor --> :yellow
      time, count = _epidemic_curve(events, State_R, min, max)
      label := "R"
      time, count
    end
  end
end

@recipe function plot(events::Events{T}) where T <: EpidemicModel
  return plot(events, minimum(events), maximum(events))
end

@recipe function plot(mcmc::MCMC{T}, mc::Int64, iter::Int64) where T <: EpidemicModel
  if !(1 <= mc <= length(mcmc.markov_chains))
    @error "Invalid Markov Chain specification"
  elseif !(1 <= iter <= mcmc.markov_chains[mc].iterations)
    @error "Invalid iteration specification"
  end
  return plot(mcmc.markov_chains[mc].events[iter])
end

@recipe function plot(network::TransmissionNetwork, pop::Population, events::Events{T},  time::Float64) where T <: EpidemicModel
  xguide --> ""
  yguide --> ""
  legend --> :topright
  xlims --> extrema(pop.risks[:x]) .+ (sum(extrema(pop.risks[:x]).*(-1,1)) .* (-0.05, 0.05))
  ylims --> extrema(pop.risks[:y]) .+ (sum(extrema(pop.risks[:y]).*(-1,1)) .* (-0.05, 0.05))
  aspect_ratio --> :equal
  ids_susceptible = _ids_by_state(events, State_S, time)
  if T in [SEIR; SEI]
    ids_exposed = _ids_by_state(events, State_E, time)
    for i in ids_exposed
      @series begin
        x, y = _plot_pathway(pop, network, i)
        seriestype := :path
        seriescolor := :black
        label := ""
        x, y
      end
    end
  end
  ids_infected = _ids_by_state(events, State_I, time)
  for i in ids_infected
    @series begin
      x, y = _plot_pathway(pop, network, i)
      seriestype := :path
      seriescolor := :black
      label := ""
      x, y
    end
  end
  if T in [SEIR; SIR]
    ids_removed = _ids_by_state(events, State_R, time)
    for i in ids_removed
      @series begin
        x, y = _plot_pathway(pop, network, i)
        seriestype := :path
        seriescolor := :black
        label := ""
        x, y
      end
    end
  end
  @series begin
    x, y = _plot_population(pop, ids_susceptible)
    seriestype := :scatter
    seriescolor --> :purple
    label := "S"
    x, y
  end
  if T in [SEIR; SEI]
    @series begin
      x, y = _plot_population(pop, ids_exposed)
      seriestype := :scatter
      seriescolor --> :lightblue4
      label := "E"
      x, y
    end
  end
  @series begin
    x, y = _plot_population(pop, ids_infected)
    seriestype := :scatter
    seriescolor --> :lightgreen
    label := "I"
    x, y
  end
  if T in [SEIR; SIR]
    @series begin
      x, y = _plot_population(pop, ids_removed)
      seriestype := :scatter
      seriescolor --> :yellow
      label := "R"
      x, y
    end
  end
end

@recipe function plot(mcmc::MCMC{T}, mc::Int64, iter::Int64, time::Float64) where T <: EpidemicModel
  if !(1 <= mc <= length(mcmc.markov_chains))
    @error "Invalid Markov Chain specification"
  elseif !(1 <= iter <= mcmc.markov_chains[mc].iterations)
    @error "Invalid iteration specification"
  end
  return plot(mcmc.markov_chains[mc].transmission_network[iter],
              mcmc.population,
              mcmc.markov_chains[mc].events[iter],
              time)
end

@recipe function plot(x::Vector{RiskParameters{T}}) where T <: EpidemicModel
  xguide --> "Iterations"
  yguide --> "Value"
  lab = [latexstring("\\Theta_{$i}") for i = 1:length(x[1])]'
  convert(Array{Float64, 2}, x)
end
