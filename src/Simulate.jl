function create_pop_db(sequence_init::Vector, covariate_init::Array)
"""
Initialize an infection database
"""
  pop_db = DataFrame(ID=0,
                     infection_times=NaN,
                     recovery_times=NaN,
                     covariate_times=NaN,
                     sequence_times=[0.],
                     infection_sources=NaN,
                     covariate_history=[[NaN NaN]],
                     sequence_history=[[sequence_init]])
  for i = 1:size(covariate_init, 1)
    push!(pop_db, DataArray(ID=i,
                            infection_times=Float64[],
                            recovery_times=Float64[],
                            covariate_times=[0.],
                            sequence_times=Float64[],
                            infection_sources=Int64[],
                            covariate_history=[covariate_init[i,:]],
                            sequence_history=Array[]))
  end
  return infection_db
end

