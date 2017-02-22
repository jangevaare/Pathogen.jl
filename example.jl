using Pathogen
using DataFrames
using Distributions
using PhyloTrees
using PhyloModels
using Plots

# Define population
x_coordinates = rand(Uniform(0, 3), 25)
y_coordinates = rand(Uniform(0, 3), 25)
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
infectivity_params = [1., 6.]
latency_params = [1/7.]
removal_params = [1/7.]

risk_params = RiskParameters(sparks_params,
                             susceptibility_params,
                             transmissibility_params,
                             infectivity_params,
                             latency_params,
                             removal_params)

# Initialize the simulation
states, rates, events, network = initialize_simulation(population,
                                                       risk_funcs,
                                                       risk_params)

# Simulate `n` events
n = 300
states, rates, events, network = simulate!(n,
                                           states,
                                           rates,
                                           events,
                                           network,
                                           population,
                                           risk_funcs,
                                           risk_params)

# Plot events
#plot(events)
#plot(population, events, network, 20.)

# Generate observations
observations = observe(events, Uniform(0., 0.5))

# Generate the associated phylogenetic tree
tree = Tree()
generatetree!(tree, events, observations, network)

# Plot the tree
#plot(tree)

# Set up a Dict for sequence data
node_data = Dict{Int64, Sequence}()

# Define a substitution model
substitution_model = JC69([1.0e-5])

# Set the root sequence
node_data[findroots(tree)[1]] = simulate(500, substitution_model)

# Simulate remaining sequences
site_rates = fill(1., 500)
simulate!(node_data, tree, substitution_model, site_rates)

# Extract observed sequences
observed_nodes = findleaves(tree)
observed_sequences = Dict{Int64, Sequence}()
for i in observed_nodes
  observed_sequences[i] = node_data[i]
end

# Set prior distributions
riskparameter_priors = RiskParameterPriors([Uniform(0., 0.001)],
                                            UnivariateDistribution[],
                                            UnivariateDistribution[],
                                            [Uniform(0., 2.), Uniform(4., 8.)],
                                            [Uniform(0., 1.)],
                                            [Uniform(0., 1.)])

substitutionmodel_priors = JC69Prior([Uniform(0., 2e-5)])

# Generate initial values for event times
# event_proposal = generate_events(observations, 10., 2., 2.)
event_proposal = events

# Initialize MCMC
phylodynamicILM_trace, phylogenetic_trace = initialize_mcmc(observations,
                                                            observed_sequences,
                                                            event_proposal,
                                                            riskparameter_priors,
                                                            risk_funcs,
                                                            substitutionmodel_priors,
                                                            population)

# Transition kernel
transition_kernel_var1 = diagm([0.000025; 0.05; 0.1; 0.025; 0.025])
transition_kernel_var2 = diagm([2.5e-7])

# Run MCMC
mcmc!(phylodynamicILM_trace,
      phylogenetic_trace,
      4999,
      transition_kernel_var1,
      transition_kernel_var2,
      1.0,
      observations,
      observed_sequences,
      riskparameter_priors,
      risk_funcs,
      substitutionmodel_priors,
      population,
      [1/4; 1/4; 1/4; 1/4])

# Tune covariance matrices
transition_kernel_var1 = cov(Array(phylodynamicILM_trace.riskparameters))/10
transition_kernel_var2 = [var([phylogenetic_trace.substitutionmodel[i].Θ[j] for i = 1:25000, j = 1])]

# Run MCMC
mcmc!(phylodynamicILM_trace,
      phylogenetic_trace,
      25000,
      transition_kernel_var1,
      diagm(transition_kernel_var2),
      0.25,
      observations,
      observed_sequences,
      riskparameter_priors,
      risk_funcs,
      substitutionmodel_priors,
      population,
      [1/4; 1/4; 1/4; 1/4])


# Tune covariance matrices
transition_kernel_var1 = transition_kernel_variance(phylodynamicILM_trace.riskparameters)/10
transition_kernel_var2 = transition_kernel_variance(phylogenetic_trace.substitutionmodel)

# Run MCMC
mcmc!(phylodynamicILM_trace,
      phylogenetic_trace,
      25000,
      transition_kernel_var1,
      diagm(transition_kernel_var2),
      0.25,
      observations,
      observed_sequences,
      riskparameter_priors,
      risk_funcs,
      substitutionmodel_priors,
      population,
      [1/4; 1/4; 1/4; 1/4])


maxiter = findlast(phylodynamicILM_trace.logposterior .== maximum(phylodynamicILM_trace.logposterior))

plot(events)
plot(phylodynamicILM_trace.events[maxiter])

plot(tree)
plot(phylogenetic_trace.tree[maxiter])

plot(population, events, network, 500.)
plot(population, phylodynamicILM_trace.events[maxiter], phylodynamicILM_trace.network[maxiter], 500.)

plot(phylodynamicILM_trace.logposterior)

plot(Array(phylodynamicILM_trace.riskparameters))
plot!([phylogenetic_trace.substitutionmodel[i].Θ[j] for i = 1:75000, j = 1])
