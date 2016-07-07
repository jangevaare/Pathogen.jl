function popplot(population::DataFrame, events::Events, time::Float64, state::Symbol)
  if !(state in [:S; :E; :I; :R])
    error("State must be specified as S, E, I or R")
  end
  x = [NaN]
  y = [NaN]
  for i in 1:length(events.susceptible)
    if findstate(events, i, time) == state
      push!(x, population[:x][i])
      push!(y, population[:y][i])
    end
  end
  return x, y
end


function pathplot(population::DataFrame, events::Events, time::Float64)
  states = findstate(events, time)
  # Count how many paths their are
  x = Array(Float64, (2,0))
  y = Array(Float64, (2,0))
  for i = 1:length(states)
    # If state is exposed, infected, or removed...
    if states[i] in [:E; :I; :R]
      # If exposure is not from external source...
      if !events.network[1][i]
        source = findfirst(events.network[2][:, i])
        x = hcat(x, population[:x][[i; source]])
        y = hcat(y, population[:y][[i; source]])
      end
    end
  end
  return x, y
end


function epiplot(events::Events, state::Symbol)
  if !(state in [:S, :E, :I, :R])
    error("State must be specified as S, E, I or R")
  end
  time=[0.]
  amount=[0]
  if state == :S
    eventorder = ind2sub((length(events.susceptible), 2), sortperm([events.susceptible events.exposed][:]))
    for i in 1:length(eventorder[1])
      if eventorder[2][i] == 1
        if isnan(events.susceptible[eventorder[1][i]])
          break
        end
        if time[end] == events.susceptible[eventorder[1][i]]
          amount[end] += 1
        else
          push!(time, events.susceptible[eventorder[1][i]])
          push!(amount, amount[end]+1)
        end
      elseif eventorder[2][i] == 2
        if isnan(events.exposed[eventorder[1][i]])
          break
        end
        if time[end] == events.exposed[eventorder[1][i]]
          amount[end] -= 1
        else
          push!(time, events.exposed[eventorder[1][i]])
          push!(amount, amount[end]-1)
        end
      end
    end
  elseif state == :E
    eventorder = ind2sub((length(events.exposed), 2), sortperm([events.exposed events.infected][:]))
    for i in 1:length(eventorder[1])
      if eventorder[2][i] == 1
        if isnan(events.exposed[eventorder[1][i]])
          break
        end
        if time[end] == events.exposed[eventorder[1][i]]
          amount[end] += 1
        else
          push!(time, events.exposed[eventorder[1][i]])
          push!(amount, amount[end]+1)
        end
      elseif eventorder[2][i] == 2
        if isnan(events.infected[eventorder[1][i]])
          break
        end
        if time[end] == events.infected[eventorder[1][i]]
          amount[end] -= 1
        else
          push!(time, events.infected[eventorder[1][i]])
          push!(amount, amount[end]-1)
        end
      end
    end
  elseif state == :I
    eventorder = ind2sub((length(events.infected), 2), sortperm([events.infected events.removed][:]))
    for i in 1:length(eventorder[1])
      if eventorder[2][i] == 1
        if isnan(events.infected[eventorder[1][i]])
          break
        end
        if time[end] == events.infected[eventorder[1][i]]
          amount[end] += 1
        else
          push!(time, events.infected[eventorder[1][i]])
          push!(amount, amount[end]+1)
        end
      elseif eventorder[2][i] == 2
        if isnan(events.removed[eventorder[1][i]])
          break
        end
        if time[end] == events.removed[eventorder[1][i]]
          amount[end] -= 1
        else
          push!(time, events.removed[eventorder[1][i]])
          push!(amount, amount[end]-1)
        end
      end
    end
  elseif state == :R
    eventorder = ind2sub((length(events.removed), 2), sortperm(events.removed))
    for i in 1:length(eventorder[1])
      if eventorder[2][i] == 1
        if isnan(events.removed[eventorder[1][i]])
          break
        end
        if time[end] == events.removed[eventorder[1][i]]
          amount[end] += 1
        else
          push!(time, events.removed[eventorder[1][i]])
          push!(amount, amount[end]+1)
        end
      end
    end
  end
  tmax = maximum([events.susceptible; events.exposed; events.infected; events.removed])
  push!(time, tmax)
  push!(amount, amount[end])
  return time, amount
end
