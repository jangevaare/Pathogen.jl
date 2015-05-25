function create_population(sequence_init::Vector, covariate_init::Array)
"""
Initialize an infection database
"""
  # infection times, infection sources, recovery times, covariate times, sequence times
  events = Array([],[],[],[],[]),

  # covariate history, sequence history
  history = Array([[ ], [ ]],[[ ], [ ]])

  return population(events, history)
end

