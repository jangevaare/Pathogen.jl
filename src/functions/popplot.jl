"""
popplot(population::DataFrame,
        events::Events,
        time::Float64,
        state::Symbol)

Return the coordinates of all individuals in a specific disease state
"""
function popplot(population::DataFrame,
                 events::Events,
                 time::Float64,
                 state::Symbol)
  if !(state in [:S; :E; :I; :R])
    throw("State must be specified as :S, :E, :I, or :R")
  end
  x = [NaN]
  y = [NaN]
  for i = 1:events.individuals
    if findstate(events, i, time) == state
      push!(x, population[:x][i])
      push!(y, population[:y][i])
    end
  end
  return x, y
end
