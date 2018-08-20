function _population_plot(population::Population,
                          ids::Vector{Int64})
  x = Float64[]
  y = Float64[]
  for i ∈ ids
    push!(x, population[:x][i])
    push!(y, population[:y][i])
  end
  return x, y
end
