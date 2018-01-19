"""
initialize_mcmc(event_obs::EventObservations,
                seq_obs::Dict{Int64, Sequence},
                events_proposal::Events,
                riskparameter_priors::RiskParameterPriors,
                riskfuncs::RiskFunctions,
                substitutionmodel_priors::SubstitutionModelPrior,
                population::DataFrame;
                conditional_network_proposals=true::Bool)

Initialize MCMC
"""
function initialize_mcmc(event_obs::EventObservations,
                         seq_obs::Dict{Int64, Sequence},
                         events_proposal::Events,
                         riskparameter_priors::RiskParameterPriors,
                         riskfuncs::RiskFunctions,
                         substitutionmodel_priors::SubstitutionModelPrior,
                         population::DataFrame;
                         conditional_network_proposals=true::Bool,
                         attempts=20::Int64)
  if attempts < 1
    error("Must make at least one initialization attempt")
  end
  initialized_trace = PathogenTrace(RiskParameters[],
                                    SubstitutionModel[],
                                    Events[],
                                    Network[],
                                    [-Inf])
  for i = 1:attempts
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
    network_proposal = propose(network_rates,
                               conditional_network_proposals = conditional_network_proposals)
    if !conditional_network_proposals
      llikelihood += loglikelihood(network_proposal, network_rates)
    end
    tree_proposal = generate_tree(events_proposal,
                                  event_obs,
                                  network_proposal)
    llikelihood += loglikelihood(tree_proposal,
                                 substitutionmodel_proposal,
                                 seq_obs)
    lposterior = lprior + llikelihood
    if lposterior > initialized_trace.logposterior[1]
      initialized_trace = PathogenTrace([riskparameter_proposal],
                                        [substitutionmodel_proposal],
                                        [events_proposal],
                                        [network_proposal],
                                        [lposterior])
    end
  end
  return initialized_trace
end
