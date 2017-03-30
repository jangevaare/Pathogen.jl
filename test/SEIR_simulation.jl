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
