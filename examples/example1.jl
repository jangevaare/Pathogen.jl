using Pathogen

JC69((2,))

init_seq = generate_sequence(10, 0.25, 0.25, 0.25, 0.25)
init_var = rand((10,2))

pop = create_population(init_seq, init_var)

pop.events
pop.history
pop.history[1][2][1]

PowerLaw = CreatePowerLaw(1., 1., 1., 1.)
PowerLaw(rand((2,10)),rand((2,10)))
