using Pathogen
using Base.Test
using DataFrames
using Distributions
using PhyloTrees

# Define population
x_coordinates = rand(Uniform(0, 10), 100)
y_coordinates = rand(Uniform(0, 10), 100)
age = rand(Poisson(40), 100)
population = DataFrame(x = x_coordinates,
                       y = y_coordinates,
                       age = age)

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

function detection_func(parameters::Vector{Float64}, population::DataFrame, k::Int64)
  return parameters[1]
end

function removal_func(parameters::Vector{Float64}, population::DataFrame, k::Int64)
  return parameters[1]
end

risk_funcs = RiskFunctions(sparks_func,
                           susceptibility_func,
                           transmissibility_func,
                           infectivity_func,
                           latency_func,
                           detection_func,
                           removal_func)

# Parametrize event rate functions
sparks_params = [0.0001]
susceptibility_params = Float64[]
transmissibility_params = Float64[]
infectivity_params = [3., 5.]
latency_params = [1/7.]
detection_params = [1.]
removal_params = [1/7.]

risk_params = RiskParameters(sparks_params,
                             susceptibility_params,
                             transmissibility_params,
                             infectivity_params,
                             latency_params,
                             detection_params,
                             removal_params)

# Initialize the simulation
index_case = 1
rates, events = initialize_simulation(population,
                                      risk_funcs,
                                      risk_params,
                                      index_case)

# Simulate `n` events
n = 300
rates, events = simulate!(n, rates, events)

# Generate the associated phylogenetic tree
tree, observed = generatetree(events)

# Simulate sequence data for each of the previously generated transmission trees
substitution_model = JC69([1.0e-5])
tree_sequences = simulate(tree, substitution_model, 1000)
