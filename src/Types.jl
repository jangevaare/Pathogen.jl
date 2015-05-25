type population
  """
  A specialized type where each individual in a population is characterized by a vector of vectors for event details and by a vector of arrays for event histories
  """
  events::Array
  history::Array
end

"""
For reference, here is the previous population database format

pop_db = DataFrame(ID=Int64,
                     infection_times=Vector[],
                     recovery_times=Vector[],
                     covariate_times=Vector[],
                     sequence_times=Vector[],
                     infection_sources=Vector[],
                     covariate_history=Array[],
                     sequence_history=Array[])
"""
