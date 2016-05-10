import Plots.scatter

function plot(population::DataFrame, events::Events, time::Float64)
  return scatter(population[:x], population[:y], color=findstate(events, 1:length(population), time))
end
