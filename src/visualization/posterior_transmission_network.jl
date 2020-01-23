function _transmission_pathway_alpha(pop::Population,
                                     network::TNPosterior,
                                     i::Int64, j::Int64)
  x = pop.risks[[i; j], :x]
  y = pop.risks[[i; j], :y]
  alpha = network.internal[j, i]
  return x, y, alpha
end

@recipe function f(tn::TNPosterior,
                   pop::Population) where T <: EpidemicModel
  xguide --> ""
  yguide --> ""
  legend --> :none
  xlims --> extrema(pop.risks[!, :x]) .+ (sum(extrema(pop.risks[!, :x]).*(-1,1)) .* (-0.05, 0.05))
  ylims --> extrema(pop.risks[!, :y]) .+ (sum(extrema(pop.risks[!, :y]).*(-1,1)) .* (-0.05, 0.05))
  aspect_ratio := :equal
  axis-->nothing
  n = pop.individuals
  for i = 1:n, j =1:n
    if tn.internal[j, i] > 0.0
      @series begin
        x, y, alpha = _transmission_pathway_alpha(pop, tn, i, j)
        seriestype --> :path
        seriescolor --> :black
        label := ""
        seriesalpha := alpha
        x, y
      end
    end
  end
end
