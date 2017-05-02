"""
pathplot(population::DataFrame,
         events::Events,
         network::Network,
         time::Float64)

Return coordinates of lines to display disease transmission network
"""
function pathplot(population::DataFrame,
                  events::Events,
                  network::Network,
                  time::Float64)
  states = [findstate(events, i, time) for i = 1:events.individuals]
  # Count how many paths their are
  x = Array(Float64, (2, 0))
  y = Array(Float64, (2, 0))
  for i = 1:length(states)
    # If state is exposed, infected, or removed...
    if states[i] in [:E; :I; :R]
      # If exposure is not from external source...
      if !network.external[i]
        source = findfirst(network.internal[:, i])
        x = hcat(x, population[:x][[i; source]])
        y = hcat(y, population[:y][[i; source]])
      end
    end
  end
  return x, y
end
