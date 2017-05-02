"""
initialize_mcmc(event_obs::EventObservations,
                seq_obs::Dict{Int64, Sequence},
                events_proposal::Events,
                riskparameter_priors::RiskParameterPriors,
                riskfuncs::RiskFunctions,
                substitutionmodel_priors::SubstitutionModelPrior,
                population::DataFrame)

Initialize MCMC
"""
function initialize_mcmc(event_obs::EventObservations,
                         seq_obs::Dict{Int64, Sequence},
                         events_proposal::Events,
                         riskparameter_priors::RiskParameterPriors,
                         riskfuncs::RiskFunctions,
                         substitutionmodel_priors::SubstitutionModelPrior,
                         population::DataFrame)
  riskparameter_proposal = propose(riskparameter_priors)
  substitutionmodel_proposal = rand(substitutionmodel_priors)
  lprior = logprior(riskparameter_priors,
                    riskparameter_proposal)
  lprior += logprior(substitutionmodel_priors,
                     substitutionmodel_proposal)
  llikelihood, network_rates = loglikelihood(riskparameter_proposal,
                                             events_proposal,
                                             riskfuncs,
                                             population)
  network_proposal = propose(network_rates)
  tree_proposal = generate_tree(events_proposal,
                                event_obs,
                                network_proposal)
  llikelihood += loglikelihood(tree_proposal,
                               substitutionmodel_proposal,
                               seq_obs)
  lposterior = lprior + llikelihood
  pathogen_trace = PathogenTrace([riskparameter_proposal],
                                 [events_proposal],
                                 [network_proposal],
                                 [lposterior])
  phylo_trace = PhyloTrace([substitutionmodel_proposal],
                           [tree_proposal],
                           [lposterior])
  return pathogen_trace, phylo_trace
end
