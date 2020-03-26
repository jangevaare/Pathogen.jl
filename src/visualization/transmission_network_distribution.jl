function _transmission_pathway_density(pop::Population,
                                       tn::TNDistribution,
                                       i::Int64, j::Int64)
  x = pop.risks[[i; j], :x]
  y = pop.risks[[i; j], :y]
  dens = tn.internal[j, i]
  return x, y, dens
end

@recipe function f(tn::TNDistribution,
                   pop::Population;
                   show_individuals=true::Bool) where M <: ILM
  xguide --> ""
  yguide --> ""
  legend --> :none
  xlims --> extrema(pop.risks[!, :x]) .+ (sum(extrema(pop.risks[!, :x]).*(-1,1)) .* (-0.05, 0.05))
  ylims --> extrema(pop.risks[!, :y]) .+ (sum(extrema(pop.risks[!, :y]).*(-1,1)) .* (-0.05, 0.05))
  aspect_ratio := :equal
  axis --> nothing
  framestyle --> :none
  linecolor --> :black
  n = pop.individuals
  for i = 1:n, j =1:n
    if tn.internal[j, i] > 0.0
      @series begin
        x, y, density = _transmission_pathway_density(pop, tn, i, j)
        seriestype --> :path
        linealpha --> density
        x, y
      end
    end
  end
  if show_individuals
    @series begin
      pop
    end
  end
end
