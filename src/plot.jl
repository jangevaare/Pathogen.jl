import Plots.plot


"""
Plot the location and disease status of a population at a specific time
"""
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
  detected_time=[0.]
  detected_amount=[0]
  removed_time=[0.]
  removed_amount=[0]
  eventorder = ind2sub((length(events.susceptible),5), sortperm([events.susceptible events.exposed events.infected events.detected events.removed][:]))
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
      if isnan(events.detected[eventorder[1][i]])
        break
      end
      if detected_time[end] == events.detected[eventorder[1][i]]
        detected_amount[end] += 1
      else
        push!(detected_time, events.detected[eventorder[1][i]])
        push!(detected_amount, detected_amount[end]+1)
      end
    elseif eventorder[2][i] == 5
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
  epidemic = DataFrame(time = [susceptible_time; exposed_time; infected_time; removed_time; detected_time],
                       count = [susceptible_amount; exposed_amount; infected_amount; removed_amount; detected_amount],
                       state = [fill("S", length(susceptible_time));
                                fill("E", length(exposed_time));
                                fill("I", length(infected_time));
                                fill("R", length(removed_time));
                                fill("D", length(detected_time))])
  return plot(epidemic,
              :time,
              :count,
              group=:state,
              xlabel="Time",
              ylabel="Amount")
end
