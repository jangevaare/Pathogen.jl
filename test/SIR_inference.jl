# Generate initial event times for MCMC based on observations
initial_events = generate_events(observations, 1., 1.)

# Set the event extents for data augmentation
event_extents = SIR_EventExtents(2.0, 2.0)

# Set prior distributions
riskparameter_priors = SIR_RiskParameterPriors([Uniform(0., 0.001)],
                                                UnivariateDistribution[],
                                                UnivariateDistribution[],
                                                [Uniform(0., 2.), Uniform(4., 8.)],
                                                [Uniform(0., 1.)])

substitutionmodel_priors = JC69Prior([Uniform(0., 2e-5)])

# Initialize MCMC
phylodynamicILM_trace, phylogenetic_trace = initialize_mcmc(observations,
                                                            observed_sequences,
                                                            initial_events,
                                                            riskparameter_priors,
                                                            risk_funcs,
                                                            substitutionmodel_priors,
                                                            population)

# Transition kernels
transition_kernel_var1 = diagm([0.0000025; 0.005; 0.01; 0.0025])
transition_kernel_var2 = diagm([2.5e-8])

# Run MCMC
mcmc!(phylodynamicILM_trace,
      phylogenetic_trace,
      1000,
      10,
      transition_kernel_var1,
      transition_kernel_var2,
      1.0,
      event_extents,
      observations,
      observed_sequences,
      riskparameter_priors,
      risk_funcs,
      substitutionmodel_priors,
      population,
      [1/4; 1/4; 1/4; 1/8; 1/8])
