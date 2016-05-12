import Plots.plot
import Plots.plot!
import Plots.ColorGradient


"""
Plot the location and disease status of a population at a specific time
"""
function plot(population::DataFrame, events::Events, time::Float64, paths=true::Bool)
  states = findstate(events, time, true)
  if paths
    # Count how many paths their are
    pathcount = 0
    # DataFrame to store exposure pathways
    pathlines = DataFrame(x = Float64[],
                          y = Float64[],
                          line = Int64[])
    for i = 1:length(states)
      # If state is exposed, infected, or removed...
      if states[i] > 1
        # If exposure is not from external source...
        if !events.network[1][i]
          pathcount += 1
          source = findfirst(events.network[2][:, i])
          append!(pathlines,
                  DataFrame(x = population[:x][[i; source]],
                            y = population[:y][[i; source]],
                            line = [i; i]))
        end
      end
    end
    # If there is at least one path to plot...
    if pathcount > 0
      plot(pathlines,
           :x,
           :y,
           group = :line,
           color = :black,
           line = (1.0, 1.0, :path),
           lab = "")
      plot!(vcat(DataFrame(x = fill(NaN, 4),
                           y = fill(NaN, 4),
                           state = ["1"; "2"; "3"; "4"]),
                 DataFrame(x = population[:x],
                           y = population[:y],
                           state = ["1"; "2"; "3"; "4"][states])),
            :x,
            :y,
            group = :state,
            line = (1, 1, :scatter),
            xlabel = "",
            ylabel = "",
            lab = ["S" "E" "I" "R"],
            seriescolor = ColorGradient(:viridis).colors[[1; 10; 20; 30]]')
    else
      plot(vcat(DataFrame(x = fill(NaN, 4),
                          y = fill(NaN, 4),
                          state = ["1", "2", "3", "4"]),
                DataFrame(x = population[:x],
                          y = population[:y],
                          state = ["1", "2", "3", "4"][states])),
           :x,
           :y,
           group = :state,
           line = (1, 1, :scatter),
           xlabel = "",
           ylabel = "",
           lab = ["S" "E" "I" "R"],
           seriescolor = ColorGradient(:viridis).colors[[1; 10; 20; 30]]')
    end
  else
    plot(vcat(DataFrame(x = fill(NaN, 4),
                        y = fill(NaN, 4),
                        state = ["1", "2", "3", "4"]),
              DataFrame(x = population[:x],
                        y = population[:y],
                        state = ["1", "2", "3", "4"][states])),
         :x,
         :y,
         group = :state,
         line = (1, 1, :scatter),
         xlabel = "",
         ylabel = "",
         lab = ["S" "E" "I" "R"],
         seriescolor = ColorGradient(:viridis).colors[[1; 10; 20; 30]]')
  end
end


"""
Epidemic plots
"""
function plot(events::Events)
  susceptible_time=[0.]
  susceptible_amount=[0]
  exposed_time=[0.]
  exposed_amount=[0]
  infected_time=[0.]
  infected_amount=[0]
  removed_time=[0.]
  removed_amount=[0]
  eventorder = ind2sub((length(events.susceptible),4), sortperm([events.susceptible events.exposed events.infected events.removed][:]))
  for i in 1:length(eventorder[1])
    if eventorder[2][i] == 1
      if isnan(events.susceptible[eventorder[1][i]])
        break
      end
      if susceptible_time[end] == events.susceptible[eventorder[1][i]]
        susceptible_amount[end] += 1
      else
        push!(susceptible_time, events.susceptible[eventorder[1][i]])
        push!(susceptible_amount, susceptible_amount[end]+1)
      end
    elseif eventorder[2][i] == 2
      if isnan(events.exposed[eventorder[1][i]])
        break
      end
      if susceptible_time[end] == events.exposed[eventorder[1][i]]
        susceptible_amount[end] -= 1
      else
        push!(susceptible_time, events.exposed[eventorder[1][i]])
        push!(susceptible_amount, susceptible_amount[end]-1)
      end
      if exposed_time[end] == events.exposed[eventorder[1][i]]
        exposed_amount[end] += 1
      else
        push!(exposed_time, events.exposed[eventorder[1][i]])
        push!(exposed_amount, exposed_amount[end]+1)
      end
    elseif eventorder[2][i] == 3
      if isnan(events.infected[eventorder[1][i]])
        break
      end
      if exposed_time[end] == events.infected[eventorder[1][i]]
        exposed_amount[end] -= 1
      else
        push!(exposed_time, events.infected[eventorder[1][i]])
        push!(exposed_amount, exposed_amount[end]-1)
      end
      if infected_time[end] == events.infected[eventorder[1][i]]
        infected_amount[end] += 1
      else
        push!(infected_time, events.infected[eventorder[1][i]])
        push!(infected_amount, infected_amount[end]+1)
      end
    elseif eventorder[2][i] == 4
      if isnan(events.removed[eventorder[1][i]])
        break
      end
      if infected_time[end] == events.removed[eventorder[1][i]]
        infected_amount[end] -= 1
      else
        push!(infected_time, events.removed[eventorder[1][i]])
        push!(infected_amount, infected_amount[end]-1)
      end
      if removed_time[end] == events.removed[eventorder[1][i]]
        removed_amount[end] += 1
      else
        push!(removed_time, events.removed[eventorder[1][i]])
        push!(removed_amount, removed_amount[end]+1)
      end
    end
  end
  tmax = maximum([susceptible_time[end]; exposed_time[end]; infected_time[end]; removed_time[end]])
  push!(susceptible_time, tmax)
  push!(exposed_time, tmax)
  push!(infected_time, tmax)
  push!(removed_time, tmax)
  push!(susceptible_amount, susceptible_amount[end])
  push!(exposed_amount, exposed_amount[end])
  push!(infected_amount, infected_amount[end])
  push!(removed_amount, removed_amount[end])

  epidemic = DataFrame(time = [susceptible_time; exposed_time; infected_time; removed_time],
                       count = [susceptible_amount; exposed_amount; infected_amount; removed_amount],
                       state = [fill("state1", length(susceptible_time));
                                fill("state2", length(exposed_time));
                                fill("state3", length(infected_time));
                                fill("state4", length(removed_time))])
  return plot(epidemic,
              :time,
              :count,
              group=:state,
              line=(0.8, 2.0, :step),
              xlabel="Time",
              ylabel="",
              lab=["S" "E" "I" "R"],
              xlim=(-1., maximum(epidemic[:time])+1),
              seriescolor = ColorGradient(:viridis).colors[[1; 10; 20; 30]]')
end
