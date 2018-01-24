# Set the event extents for generation of the initial event times
initial_event_extents = SI_EventExtents(1.)

# Set prior distributions
riskparameter_priors = SI_RiskParameterPriors([Uniform(0., 0.001)],
                                              UnivariateDistribution[],
                                              UnivariateDistribution[],
                                              [Uniform(0., 2.), Uniform(4., 8.)])

substitutionmodel_priors = JC69Prior([Uniform(0., 2e-5)])

# Initialize MCMC
phylodynamicILM_trace = initialize_mcmc(observations,
                                        initial_event_extents,
                                        observed_sequences,
                                        riskparameter_priors,
                                        risk_funcs,
                                        substitutionmodel_priors,
                                        population)

# Transition kernels
transition_kernel = [0.0000025; 0.005; 0.01; 2.5e-8]

# Set the event extents for data augmentation
event_extents = SI_EventExtents(2.)

# Run MCMC
mcmc!(phylodynamicILM_trace,
      100,
      2,
      transition_kernel,
      1.0,
      event_extents,
      observations,
      observed_sequences,
      riskparameter_priors,
      risk_funcs,
      substitutionmodel_priors,
      population)
