using Pathogen

init_seq = generate_sequence(10, 0.25, 0.25, 0.25, 0.25)
init_var = rand((2,10))

pop = create_population(init_seq, init_var)

powerlaw = create_powerlaw(1., 1., 1., 1.)
latency = create_constantrate(3.)
recovery = create_constantrate(5.)

ratearray = create_ratearray(pop, powerlaw, jc69((0.5,)))
ratearray.rates
ratearray.events

while pop.timeline[1][end] < 10.
  onestep!(ratearray, pop, powerlaw, latency, recovery, jc69((2,)))
end

