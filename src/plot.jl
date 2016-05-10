import Plots.plot

function plot(population::DataFrame, events::Events, time::Float64)
  return plot(vcat(DataFrame(x=fill(NaN, 4), y=fill(NaN, 4), state=["S", "E", "I", "R"]),
                   DataFrame(x = population[:x], y = population[:y], state=findstate(events, time))),
                   :x,
                   :y,
                   group=:state,
                   line=(1, 1, :scatter),
                   xlabel="",
                   ylabel="",
                   title="Time = $time")
end
