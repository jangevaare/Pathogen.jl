using Pathogen

init_seq = GenerateSequence(10, 0.25, 0.25, 0.25, 0.25)
init_var = rand((2,10))

pop = CreatePopulation(init_seq, init_var)

pop.events
pop.history

PowerLaw = CreatePowerLaw(1., 1., 1., 1.)
Latency = CreateConstantRate(1.)
Recovery = CreateConstantRate(1.)

PowerLaw(pop, 2, 4)

CreateRateArray(pop, PowerLaw, Latency, Recovery, JC69((2,)))
