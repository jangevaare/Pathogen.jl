# Define event rate functions
function sparks_func(parameters::Vector{Float64}, population::DataFrame, i::Int64)
  return parameters[1]
end

function susceptibility_func(parameters::Vector{Float64}, population::DataFrame, i::Int64)
  return 1.
end

function transmissibility_func(parameters::Vector{Float64}, population::DataFrame, k::Int64)
  return 1.
end

function infectivity_func(parameters::Vector{Float64}, population::DataFrame, i::Int64, k::Int64)
  α = parameters[1]
  β = parameters[2]
  d = sqrt((population[:x][i] - population[:x][k])^2 + (population[:y][i] - population[:y][k])^2)
  return α * d^-β
end

function latency_func(parameters::Vector{Float64}, population::DataFrame, j::Int64)
  return parameters[1]
end

function removal_func(parameters::Vector{Float64}, population::DataFrame, k::Int64)
  return parameters[1]
end


# Parametrize event rate functions
sparks_params = [0.0001]
susceptibility_params = Float64[]
transmissibility_params = Float64[]
infectivity_params = [1., 6.]
latency_params = [1/15.]
removal_params = [1/3.]
