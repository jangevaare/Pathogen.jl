function create_population(init_seq::Vector, init_var::Array)
"""
Initialize an infection database. `init_seq` is assigned to the "external" infection source. Each row of the `init_array` is assigned to a individual
"""
  # identifier, infection times, infection sources, recovery times, covariate times, sequence times
  events = Array[Array[[0], [NaN],  [NaN],  [NaN],  [NaN], [0]],
                 Array[[1], [],     [],     [],     [0],   []]]

  # covariate history, sequence history
  history = Array[Array[[[fill(NaN, length(init_var[1,:]))]],     [init_seq]],
                  Array[[[init_var[1,:]]],[]]]

  # push individuals to these arrays.
  for r = 2:shape(init_var,1)
    push!(events, Array[[r], [],     [],     [],     [0],   []])
    push!(history, Array[[[init_var[r,:]]],[]])
  end

  # save as a population object type
  return population(events, history)
end

