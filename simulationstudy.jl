using Pathogen
using PhyloTrees
using PhyloModels
using DataFrames
using JLD
using Distributions

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

# 10 replicates
for i = 1:10
  # Define population
  x_coordinates = rand(Uniform(0, 5), 25)
  y_coordinates = rand(Uniform(0, 5), 25)
  population = DataFrame(x = x_coordinates,
                         y = y_coordinates)

  # 4 scenarios
  for j = 1:4
    println("Starting replicate $i of scenario $j")
    # Parametrize event rate functions
    sparks_params = [0.0001]
    susceptibility_params = Float64[]
    transmissibility_params = Float64[]
    if j == 1
      infectivity_params = [2., 7.]
      latency_params = [1/20.]
    elseif j == 2
      infectivity_params = [1., 5.]
      latency_params = [1/20.]
    elseif j == 3
      infectivity_params = [2., 7.]
      latency_params = [1/5.]
    elseif j == 4
      infectivity_params = [1., 5.]
      latency_params = [1/5.]
    end
    removal_params = [1/3.]

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

  # Generate observations
  observations = observe(events, Uniform(0., 0.5))

  # Generate the associated phylogenetic tree
  tree = generate_tree(events, observations, network)

  # Set up a Dict for sequence data
  node_data = Dict{Int64, Sequence}()

  # Define a substitution model
  substitution_model = JC69([1.0e-5])

  # Set the root sequence
  node_data[findroots(tree)[1]] = simulate(1000, substitution_model)

  # Simulate remaining sequences
  site_rates = fill(1., 1000)
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
                                              [Uniform(0., 3.), Uniform(3., 9.)],
                                              [Uniform(0., 0.25)],
                                              [Uniform(0., 1.)])

  substitutionmodel_priors = JC69Prior([Uniform(0., 2e-5)])

  # Generate initial values for event times
  # event_proposal = generate_events(observations, 21., 2., 2.)
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
  transition_kernel_var1 = diagm([0.0000025; 0.005; 0.01; 0.000625; 0.0025])
  transition_kernel_var2 = diagm([2.5e-8])

  # Run MCMC
  mcmc!(phylodynamicILM_trace,
        phylogenetic_trace,
        50000,
        10,
        transition_kernel_var1,
        transition_kernel_var2,
        1.0,
        Inf,
        2.,
        2.,
        observations,
        observed_sequences,
        riskparameter_priors,
        risk_funcs,
        substitutionmodel_priors,
        population,
        [1/4; 1/4; 1/4; 1/8; 1/8])

  # Tune covariance matrices
  transition_kernel_var1 = cov(Array(phylodynamicILM_trace.riskparameters))/10
  transition_kernel_var2 = [var([phylogenetic_trace.substitutionmodel[i].Θ[j] for i = 1:length(phylogenetic_trace), j = 1])]

  starttime = now()
  # Run MCMC
  mcmc!(phylodynamicILM_trace,
        phylogenetic_trace,
        10000000,
        250,
        transition_kernel_var1,
        diagm(transition_kernel_var2),
        1.0,
        Inf,
        2.,
        2.,
        observations,
        observed_sequences,
        riskparameter_priors,
        risk_funcs,
        substitutionmodel_priors,
        population,
        [1/4; 1/4; 1/4; 1/8; 1/8])
elapsedtime = now() - starttime

save("/home/jangevaare/Documents/Pathogen/simulation_replicate$i _scenario$j.jld",
     "phylodynamicILM_trace", phylodynamicILM_trace[5002:45001],
     "phylogenetic_trace", phylogenetic_trace[5002:45001],
     "events", events,
     "observations", observations,
     "observed_sequences", observed_sequences,
     "network", network,
     "population", population,
     "elapsedtime", elapsedtime)
  end
end
