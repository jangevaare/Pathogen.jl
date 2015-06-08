using Pathogen

init_seq = generate_sequence(10, 0.25, 0.25, 0.25, 0.25)
init_var = rand((2,10))

pop = create_population(init_seq, init_var)

powerlaw = create_powerlaw(1., 1., 1., 1.)
latency = create_constantrate(1.)
recovery = create_constantrate(1.)

ratearray = create_ratearray(pop, powerlaw, jc69((2,)))



ratearray.rates
ratearray.events


powerlaw(pop, 2, 1)
