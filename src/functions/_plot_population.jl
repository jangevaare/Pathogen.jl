function _plot_population(pop::Population,
                          ids::Vector{Int64})
  x = Float64[]
  y = Float64[]
  for i in ids
    push!(x, pop.risks[:x][i])
    push!(y, pop.risks[:y][i])
  end
  return x, y
end
