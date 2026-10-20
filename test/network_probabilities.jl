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
