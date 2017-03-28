"""
epiplot(events::Events,
        state::Symbol)

Return counts of individuals in a specified disease state and the associated
time point
"""
function epiplot(events::Events,
                 state::Symbol)
  if !(state in [:S, :E, :I, :R])
    throw("State must be specified as :S, :E, :I or :R")
  end
  time=[0.]
  amount=[0]
  if state == :S
    # All individuals start as susceptible
    amount[end] += events.individuals
    eventorder = sortperm(events.exposed)
    for i in 1:length(eventorder)
      if isnan(events.exposed[eventorder[i]])
        break
      end
      if time[end] == events.exposed[eventorder[i]]
        amount[end] -= 1
      else
        push!(time, events.exposed[eventorder[i]])
        push!(amount, amount[end]-1)
      end
    end
  elseif state == :E
    eventorder = ind2sub((events.individuals, 2), sortperm([events.exposed events.infected][:]))
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
    eventorder = ind2sub((events.individuals, 2), sortperm([events.infected events.removed][:]))
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
    eventorder = ind2sub((events.individuals, 2), sortperm(events.removed))
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
  tmax = maximum([events.exposed; events.infected; events.removed])
  push!(time, tmax)
  push!(amount, amount[end])
  return time, amount
end