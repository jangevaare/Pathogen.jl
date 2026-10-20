@testset "Truncated event proposal densities" begin
  obs = EventObservations{SI}([1.0])
  extents = EventExtents{SI}(1.0)
  current, proposed = Events{SI}([0.05]), Events{SI}([0.5])
  σ = 0.3
  forward = Pathogen._event_logproposal(proposed, current, [1], σ, extents, obs)
  reverse = Pathogen._event_logproposal(current, proposed, [1], σ, extents, obs)
  normal = Normal()
  mass(t) = cdf(normal, (1.0-t)/σ) - cdf(normal, -t/σ)
  @test reverse - forward ≈ log(mass(0.05)) - log(mass(0.5))
  @test reverse - forward < 0.0
  # Uniform targets must satisfy detailed balance even close to a bound.
  @test exp(forward) * min(1.0, exp(reverse-forward)) ≈
        exp(reverse) * min(1.0, exp(forward-reverse))
  @test Pathogen._event_logproposal(Events{SI}([1.1]), current, [1], σ, extents, obs) == -Inf
  unchanged, log_forward = Pathogen._propose_events(MersenneTwister(1), current, [1], 0.0, extents, obs)
  @test unchanged.infection == current.infection
  @test unchanged.infection !== current.infection
  @test log_forward == Pathogen._event_logproposal(current, unchanged, [1], 0.0, extents, obs) == 0.0
  @test Pathogen._event_logproposal(proposed, current, [1], 0.0, extents, obs) == -Inf
end

@testset "Batch proposal densities follow changing bounds" begin
  obs = EventObservations{SEI}([1.0])
  extents = EventExtents{SEI}(1.0, 1.0)
  current = Events{SEI}([0.2], [0.8])
  proposed = Events{SEI}([0.6], [0.9])
  σ = 0.3
  # Exposure first: infection is bounded below by the newly proposed exposure.
  expected_forward = logpdf(truncated(Normal(0.2, σ), -0.2, 0.8), 0.6) +
                     logpdf(truncated(Normal(0.8, σ), 0.6, 1.0), 0.9)
  expected_reverse = logpdf(truncated(Normal(0.6, σ), -0.1, 0.9), 0.2) +
                     logpdf(truncated(Normal(0.9, σ), 0.2, 1.0), 0.8)
  @test Pathogen._event_logproposal(proposed, current, [1, 2], σ, extents, obs) ≈ expected_forward
  @test Pathogen._event_logproposal(current, proposed, [1, 2], σ, extents, obs) ≈ expected_reverse
  sampled, log_forward = Pathogen._propose_events(MersenneTwister(2), current, [1, 2], σ, extents, obs)
  @test log_forward ≈ Pathogen._event_logproposal(sampled, current, [1, 2], σ, extents, obs)
  @test current.exposure == [0.2] && current.infection == [0.8]
  # A valid forward path can have zero reverse support in the same index order.
  x, y = Events{SEI}([0.2], [0.4]), Events{SEI}([0.6], [0.8])
  @test isfinite(Pathogen._event_logproposal(y, x, [2, 1], σ, extents, obs))
  @test Pathogen._event_logproposal(x, y, [2, 1], σ, extents, obs) == -Inf

  # In coupled inference the proposed network determines the reverse bounds.
  obs = EventObservations{SIR}([0.1, 1.0], [1.0, 2.0])
  extents = EventExtents{SIR}(1.0, 2.0)
  network = TransmissionNetwork(2)
  network.external[1] = true
  Pathogen.update!(network, EndogenousTransmission(2, 1))
  target = Events{SIR}([0.1, 0.2], [0.4, 1.5])
  source = Events{SIR}([0.1, 0.6], [0.8, 1.5])
  @test Pathogen._event_logproposal(target, source, [2], σ, extents, obs, network) ≈
        logpdf(truncated(Normal(0.6, σ), 0.1, 0.8), 0.2)
end

function uniform_time_chain(phylo; coupled=false, draws=6000)
  pop = Population(DataFrame(x=[1.0]), zeros(1, 1))
  rf = RiskFunctions{SI}(perf_constant, perf_one, perf_distance_pair, perf_one)
  rp = RiskParameters{SI}([1.5], Float64[], Float64[], Float64[])
  priors = RiskPriors{SI}(UnivariateDistribution[Uniform(1.0, 2.0)],
    UnivariateDistribution[], UnivariateDistribution[], UnivariateDistribution[])
  extents = EventExtents{SI}(1.0)
  events = Events{SI}([0.5])
  network = TransmissionNetwork(1)
  network.external[1] = true
  if phylo
    obs = EventObservations{SI}([1.0], Union{Nothing, RNASeq}[RNASeq("AC")])
    sm = JC69(0.05)
    mcmc = MCMC(obs, extents, pop, rf, priors, JC69abs, UnivariateDistribution[Uniform(0.01, 0.1)])
    tree, nodes = generate(Tree, events, obs, network)
    log_posterior = Pathogen.PhyloModels.loglikelihood(tree, sm, Dict(nodes[1] => obs.seq[1])) + log(1/0.09)
    chain = Pathogen.MarkovChain{SI, PhyloILM}(events, network, rp, sm, log_posterior)
  else
    obs = EventObservations{SI}([1.0])
    mcmc = MCMC(obs, extents, pop, rf, priors)
    chain = Pathogen.MarkovChain{SI, TNILM}(events, network, rp, 0.0)
  end
  rng = MersenneTwister(28)
  cache, work = Pathogen.TransmissionRateCache(1), Pathogen.LikelihoodWorkspace(SI, 1)
  for i in 1:draws
    if phylo
      Pathogen.update!(rng, chain, mcmc, fill(0.01, 1, 1), zeros(0, 0), 0.35;
        DA_couple=coupled, transmission_rate_cache=cache, likelihood_workspace=work)
    else
      Pathogen.update!(rng, chain, mcmc, fill(0.01, 1, 1), 0.35;
        transmission_rate_cache=cache, likelihood_workspace=work)
    end
  end
  return [e.infection[1] for e in chain.events[1001:end]]
end

@testset "Event augmentation preserves a uniform time target" begin
  # For one SI individual, conditioning on the first infection leaves a uniform
  # time target on [0,1]. One observed sequence also has a constant marginal
  # likelihood under JC69, so the same check exercises both PhyloILM modes.
  for (phylo, coupled) in ((false, false), (true, false), (true, true))
    times = uniform_time_chain(phylo; coupled=coupled)
    @test all(t -> 0.0 <= t <= 1.0, times)
    @test abs(count(t -> t < 0.1 || t > 0.9, times)/length(times) - 0.2) < 0.035
    @test abs(mean(times) - 0.5) < 0.03
  end
end
