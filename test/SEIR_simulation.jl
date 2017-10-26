risk_funcs = SEIR_RiskFunctions(sparks_func,
                                susceptibility_func,
                                transmissibility_func,
                                infectivity_func,
                                latency_func,
                                removal_func)

risk_params = SEIR_RiskParameters(sparks_params,
                                  susceptibility_params,
                                  transmissibility_params,
                                  infectivity_params,
                                  latency_params,
                                  removal_params)

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
observations = observe(events, Uniform(0., 0.5), Uniform(0., 0.5))

# Generate the associated phylogenetic tree
trees, tree_ids, leaf_ids = generate_tree(events, observations, network)

# Set up a vector of Dicts for sequence data
tree_data = fill(Dict{Int64, Sequence}(), length(trees))

# Define a substitution model
substitution_model = JC69([1.0e-4])
SNPs = 50
site_rates = fill(1., SNPs)

@simd for i = 1:length(trees)
    # Set the root sequences
    tree_data[i][1] = simulate(SNPs, substitution_model)
    # Simulate remaining sequences
    simulate!(tree_data[i], trees[i], substitution_model, site_rates)
end

# Set up a vector of Dicts for observed sequences
observed_sequences = Dict{Int64, Sequence}()

# Extract observed sequences
for i in 1:observations.individuals
    if !isnan(observations.infected[i])
        observed_sequences[i] = tree_data[tree_ids[i]][leaf_ids[i]]
    end
end
