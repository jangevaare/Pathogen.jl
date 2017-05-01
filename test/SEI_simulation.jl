risk_funcs = SEI_RiskFunctions(sparks_func,
                               susceptibility_func,
                               transmissibility_func,
                               infectivity_func,
                               latency_func)

risk_params = SEI_RiskParameters(sparks_params,
                                 susceptibility_params,
                                 transmissibility_params,
                                 infectivity_params,
                                 latency_params)

states, rates, events, network = initialize_simulation(population,
                                                       risk_funcs,
                                                       risk_params)

states, rates, events, network = simulate!(50,
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
