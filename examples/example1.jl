using Pathogen, Gadfly, Distributions, BioSeq

generate_seq(200, 0.25, 0.25, 0.25, 0.25)
generate_2bitseq(200, 0.25, 0.25, 0.25, 0.25)

init_var = rand(Uniform(0,25), (2,200))

recovery = create_constantrate(1/5.)
substitution = jc69((0.1,))

ratearray = create_ratearray(pop, powerlaw, substitution)
ratearray.rates
ratearray.events

@time while length(pop.timeline[1]) < 30000.
  onestep!(ratearray, pop, powerlaw, latency, recovery, substitution)
end

images = 1000
# Plot it
for time = 1:images
  states, routes = plotdata(pop, (time*pop.timeline[1][end])/images)
  p1 = plot(layer(states, x="x", y="y", color="state", Geom.point),
            layer(routes, x="x", y="y", group="age", Geom.polygon),
            Theme(panel_opacity=1., panel_fill=color("white"), default_color=color("black"), background_color=color("white")))
  filenumber = time/images
  filenumber = prod(split("$filenumber", ".", 2))
  filenumber *= prod(fill("0", 5-length(filenumber)))
  draw(PNG(homedir()"/Desktop/plots/infection_$filenumber.png", 15cm, 10cm), p1)
end
