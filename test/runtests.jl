using Pathogen
using DataFrames
using Distributions
using PhyloTrees
using Plots

# Define population
x_coordinates = rand(Uniform(0, 9), 100)
y_coordinates = rand(Uniform(0, 9), 100)
population = DataFrame(x = x_coordinates,
                       y = y_coordinates)

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

risk_funcs = RiskFunctions(sparks_func,
                           susceptibility_func,
                           transmissibility_func,
                           infectivity_func,
                           latency_func,
                           removal_func)

# Parametrize event rate functions
sparks_params = [0.0001]
susceptibility_params = Float64[]
transmissibility_params = Float64[]
infectivity_params = [3., 5.]
latency_params = [1/7.]
removal_params = [1/7.]

risk_params = RiskParameters(sparks_params,
                             susceptibility_params,
                             transmissibility_params,
                             infectivity_params,
                             latency_params,
                             removal_params)

# Initialize the simulation
rates, events, network = initialize_simulation(population,
                                               risk_funcs,
                                               risk_params)

# Simulate `n` events
n = 250
rates, events, network = simulate!(n, rates, events, network)

# Plot events
plot(events)
plot(population, events, network, 20.)

# Generate observations
observations = observe(events, Uniform(0., 0.5))

# Generate the associated phylogenetic tree
tree = generatetree(events, observations, network)

# Simulate sequence data for each of the previously generated transmission trees
substitution_model = JC69([1.0e-5])
observed_sequences = simulate(tree, substitution_model, 200)[findleaves(tree)]

plot(tree)

# Inference
riskparameter_priors = RiskParameterPriors([Uniform(0., 0.001)],
                                            UnivariateDistribution[],
                                            UnivariateDistribution[],
                                            [Uniform(1., 5.), Uniform(3., 7.)],
                                            [Uniform(0., 1.)],
                                            [Uniform(0., 1.)])

# Generate prior distributions for event times
event_priors = generate_eventpriors(observations, 7., 3., 3.)

# Substitution model priors
substitutionmodel_priors = JC69Prior([Uniform(5e-6, 2e-5)])

# Transition kernel
transition_kernel_var1 = transition_kernel_variance(riskparameter_priors)
transition_kernel_var2 = transition_kernel_variance(substitutionmodel_priors)

# Run MCMC
phylodynamicILM_trace, phylogenetic_trace = mcmc(1000,
                                                 transition_kernel_var1,
                                                 transition_kernel_var2,
                                                 observations,
                                                 observed_sequences,
                                                 event_priors,
                                                 riskparameter_priors,
                                                 risk_funcs,
                                                 substitutionmodel_priors,
                                                 population)
