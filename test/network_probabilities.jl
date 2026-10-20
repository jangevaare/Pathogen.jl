@testset "Conditional transmission network probabilities" begin
  tr = Pathogen.TransmissionRates([0.0, 2.0, 1.0], [0.0 2.0 3.0; 0.0 0.0 1.0; 0.0 0.0 0.0])
  states = [State_I, State_S, State_S]
  network = TransmissionNetwork(states)
  Pathogen.update!(network, EndogenousTransmission(2, 1))
  Pathogen.update!(network, EndogenousTransmission(3, 1))
  score(net, prior=nothing) = Pathogen._network_logprobabilities(net, tr, prior, [2, 3], states)
  @test score(network) ≈ [0.0, log(1/2), log(3/5)]
  prior = TNDistribution([1.0, 3.0, 2.0], [0.0 1.0 2.0; 0.0 0.0 4.0; 0.0 0.0 0.0])
  @test score(network, prior) ≈ [0.0, log(1/4), log(1/2)]
  # Enumerate every possible parent assignment, including external sources.
  mass = 0.0
  weighted_mass = 0.0
  for p2 in (0, 1), p3 in (0, 1, 2)
    net = TransmissionNetwork(states)
    Pathogen.update!(net, p2 == 0 ? ExogenousTransmission(2) : EndogenousTransmission(2, p2))
    Pathogen.update!(net, p3 == 0 ? ExogenousTransmission(3) : EndogenousTransmission(3, p3))
    mass += exp(sum(score(net)))
    weighted_mass += exp(sum(score(net, prior)))
  end
  @test mass ≈ 1.0
  @test weighted_mass ≈ 1.0
  missing = TransmissionNetwork(states)
  @test sum(score(missing)) == -Inf
  multiple = copy(network)
  multiple.external[2] = true
  @test sum(score(multiple)) == -Inf
  zero = TNDistribution(zeros(3), zeros(3, 3))
  @test sum(score(network, zero)) == -Inf
  @test !Pathogen._has_network_support(tr, zero, [2, 3])
  invalid_initial = copy(network)
  invalid_initial.external[1] = false
  @test sum(score(invalid_initial)) == -Inf
  # The initial infection is conditioned upon, not sampled as an extra source.
  @test score(network)[1] == 0.0
  survivor = copy(network)
  @test Pathogen._network_logprobabilities(survivor, tr, nothing, [2], states)[3] == -Inf
  proposed, current = [0.0, log(1/4), log(1/2)], [0.0, log(1/2), log(3/5)]
  @test Pathogen._retained_network_logratio(proposed, current, [2]) ≈ log(1/2) - log(3/5)
  @test Pathogen._retained_network_logratio(proposed, current, [2, 3]) == 0.0
end

function joint_network_fixture(phylo)
  pop = Population(DataFrame(x=ones(3)), [0.0 1.0 2.0; 0.0 0.0 4.0; 0.0 0.5 0.0])
  rf = RiskFunctions{SI}((p, pop, i) -> 0.0, perf_one, perf_distance_pair, perf_one)
  rp = RiskParameters{SI}((Float64[] for _ in 1:4)...)
  priors = RiskPriors{SI}((UnivariateDistribution[] for _ in 1:4)...)
  states = [State_I, State_S, State_S]
  prior = TNDistribution(zeros(3), [0.0 1.0 2.0; 0.0 0.0 1.0; 0.0 3.0 0.0])
  times = [-Inf, 1.0, 1.0]
  obs = phylo ? EventObservations{SI}(times, Union{Nothing, RNASeq}[RNASeq(""), RNASeq(""), RNASeq("")]) :
                EventObservations{SI}(times)
  extents = EventExtents{SI}(1.0)
  mcmc = phylo ? MCMC(obs, extents, pop, states, rf, priors, JC69abs,
                      UnivariateDistribution[Uniform(0.01, 0.1)]; tnprior=prior) :
                 MCMC(obs, extents, pop, states, rf, priors; tnprior=prior)
  return mcmc, rp
end

# Recompute the joint score independently: replay hazards, then look up each
# recorded parent's rate and normalize the prior-weighted source column.
function reference_joint_score(mc, mcmc, index)
  events, rp = mc.events[index], mc.risk_parameters[index]
  ll, rates, transmissions = Pathogen.loglikelihood(rp, mcmc.risk_functions, events,
    mcmc.population, mcmc.starting_states)
  base = ll + Pathogen.logprior(rp, mcmc.risk_priors)
  if mc.substitution_model !== nothing
    sm = mc.substitution_model[index]
    tree, nodes = generate(Tree, events, mcmc.event_observations, mc.transmission_network[index])
    data = Dict(nodes[i] => mcmc.event_observations.seq[i] for i in eachindex(nodes)
                if nodes[i] !== nothing && mcmc.event_observations.seq[i] !== nothing)
    base += Pathogen.logprior(sm, mcmc.substitution_model_prior) +
            Pathogen.PhyloModels.loglikelihood(tree, sm, data)
  end
  net, prior = mc.transmission_network[index], mcmc.transmission_network_prior
  source_score = 0.0
  for id in transmissions
    weights = [rates.external[id]; rates.internal[:, id]]
    if prior !== nothing
      weights .*= [prior.external[id]; prior.internal[:, id]]
    end
    parent = net.external[id] ? 1 : findfirst(net.internal[:, id]) + 1
    source_score += log(weights[parent] / sum(weights))
  end
  return base + source_score
end

@testset "Joint network target across inference modes" begin
  # Two finite infection times lie in [0,1]. The density when 2 precedes 3
  # is 2*exp(-6*(t3-t2)); in the opposite order it is exp(-1.5*(t2-t3)).
  # Integrating each triangle gives the exact marginal ordering probability.
  integral(k) = (k - 1 + exp(-k)) / k^2
  p2first = 2integral(6) / (2integral(6) + integral(1.5))
  for mode in (:collapsed, :conditioned, :separate, :coupled)
    phylo = mode in (:separate, :coupled)
    mcmc, rp = joint_network_fixture(phylo)
    rng = MersenneTwister(493)
    mc = Pathogen._initialization_attempt(rng, mcmc, -Inf)
    @test mc !== nothing
    @test mc.log_posterior[1] ≈ reference_joint_score(mc, mcmc, 1)
    work = Pathogen.LikelihoodWorkspace(SI, 3)
    cache = Pathogen.TransmissionRateCache(3)
    for i in 1:16000
      if phylo
        Pathogen.update!(rng, mc, mcmc, zeros(0, 0), zeros(0, 0), 0.4;
          event_batches=2, DA_couple=mode == :coupled,
          likelihood_workspace=work, transmission_rate_cache=cache)
      else
        # Also exercise switching between the collapsed and fixed-network modes.
        conditioned = mode == :conditioned && (i <= 8000 || isodd(i))
        Pathogen.update!(rng, mc, mcmc, zeros(0, 0), 0.4;
          event_batches=2, condition_on_network=conditioned,
          likelihood_workspace=work, transmission_rate_cache=cache)
      end
    end
    indices = 2002:length(mc.events)
    order_probability = mean(mc.events[i].infection[2] < mc.events[i].infection[3] for i in indices)
    @test abs(order_probability - p2first) < 0.035
    # Prior-weighted source probabilities are 1/2 for 2->3 when 2 is first,
    # and 3/5 for 3->2 when 3 is first. Empty sequences leave the epidemic target unchanged.
    @test abs(mean(mc.transmission_network[i].internal[2, 3] for i in indices) - p2first/2) < 0.035
    @test abs(mean(mc.transmission_network[i].internal[3, 2] for i in indices) - (1-p2first)*3/5) < 0.035
    for i in (1, 4000, 8000, 16001)
      @test mc.log_posterior[i] ≈ reference_joint_score(mc, mcmc, i)
    end
    @test mc.log_posterior[end] ≈ mc.log_base_posterior + sum(mc.network_logprobabilities)
  end
end

@testset "Parameter moves retain probabilities of fixed parents" begin
  mcmc, _ = joint_network_fixture(true)
  mcmc.transmission_network_prior = nothing
  # Only source 2's infectivity varies. Parent 3 stays source 1, so its
  # epidemic hazard (2+4θ) cancels its source probability 2/(2+4θ).
  mcmc.risk_functions = RiskFunctions{SI}((p, pop, i) -> 0.0, perf_one,
    (p, pop, i, j) -> pop.distances[j, i] * (j == 2 ? p[1] : 1.0), perf_one)
  mcmc.risk_priors = RiskPriors{SI}(UnivariateDistribution[], UnivariateDistribution[],
    UnivariateDistribution[Uniform(0.1, 4.0)], UnivariateDistribution[])
  rp = RiskParameters{SI}(Float64[], Float64[], [0.5], Float64[])
  events = Events{SI}([-Inf, 0.25, 0.75])
  network = TransmissionNetwork(mcmc.starting_states)
  Pathogen.update!(network, EndogenousTransmission(2, 1))
  Pathogen.update!(network, EndogenousTransmission(3, 1))
  sm = JC69(0.05)
  base = Pathogen.loglikelihood(rp, mcmc.risk_functions, events, mcmc.population,
                                mcmc.starting_states)[1] +
         Pathogen.logprior(rp, mcmc.risk_priors) + Pathogen.logprior(sm, mcmc.substitution_model_prior)
  mc = Pathogen.MarkovChain{SI, PhyloILM}(events, network, rp, sm, base)
  rng = MersenneTwister(918)
  cache, work = Pathogen.TransmissionRateCache(3), Pathogen.LikelihoodWorkspace(SI, 3)
  for i in 1:14000
    Pathogen.update!(rng, mc, mcmc, fill(0.3^2, 1, 1), zeros(0, 0), 0.0;
      event_batches=0, DA_couple=true, transmission_rate_cache=cache, likelihood_workspace=work)
  end
  # A truncated exponential density exp(-2θ) is the exact conditional target.
  expected = (exp(-0.2) - exp(-2.0)) / (exp(-0.2) - exp(-8.0))
  @test abs(mean(p.infectivity[1] < 1.0 for p in mc.risk_parameters[2002:end]) - expected) < 0.035
  @test mc.log_posterior[end] ≈ reference_joint_score(mc, mcmc, length(mc.events))
  @test all(net -> net.internal == network.internal, mc.transmission_network)
end

@testset "Network support during initialization" begin
  mcmc, _ = joint_network_fixture(false)
  mcmc.transmission_network_prior = TNDistribution(zeros(3), zeros(3, 3))
  @test Pathogen._initialization_attempt(MersenneTwister(3), mcmc, -Inf) === nothing
end
