"""
Example 1
SEIR simulation, visualization, and inference
Justin Angevaare
"""

using Pathogen, Gadfly, DataFrames, Distributions

init_seq = create_seq(200, 0.25, 0.25, 0.25, 0.25)
init_var = rand(Uniform(0,25), (2,200))

pop = create_population(init_seq, init_var)

powerlaw = create_powerlaw(4., 5., 0.002)
latency = create_constantrate(1/3.)
recovery = create_constantrate(1/5.)
substitution = jc69((0.1,))

ratearray = create_ratearray(pop, powerlaw, substitution)

ratearray.rates
ratearray.events

@time while length(pop.timeline[1]) < 30000.
  onestep!(ratearray, pop, powerlaw, latency, recovery, substitution)
end

# Simulation visualization
# images = 1000
# for time = 1:images
#   states, routes = plotdata(pop, (time*pop.timeline[1][end])/images)
#   p1 = plot(layer(states, x="x", y="y", color="state", Geom.point),
#             layer(routes, x="x", y="y", group="age", Geom.polygon),
#             Theme(panel_opacity=1., panel_fill=color("white"), default_color=color("black"), background_color=color("white")))
#   filenumber = time/images
#   filenumber = prod(split("$filenumber", ".", 2))
#   filenumber *= prod(fill("0", 5-length(filenumber)))
#   draw(PNG(homedir()"/Desktop/plots/infection_$filenumber.png", 15cm, 10cm), p1)
# end

# Inference
actual, obs = SEIR_surveilance(pop, 2.)

"""
α, β: powerlaw exposure kernel parameters
η: external pressure rate
ρ: infectivity rate (1/mean latent period)
γ: recovery rate (1/mean infectious period)
ν: detection rate (1/mean detection lag)
"""

#priors = SEIR_priors(Gamma(8.,2.), Gamma(10.,2.), Gamma(0.5,10.), Gamma(1.,3.), Gamma(1.,5.), Gamma(2.,1.))
priors = SEIR_priors(Uniform(0,10), Uniform(0,10), Uniform(0,0.1), Uniform(0,1), Uniform(0,1), Uniform(0,5))

trace = SEIR_initialize(priors, obs)

@time SEIR_MCMC(100000, diagm([0.5, 0.5, 0.5, 0.5, 0.5, 0.5]), trace, priors, obs)

for i = 1:16
  opt_cov = cov([trace.α trace.β trace.ρ trace.γ trace.η trace.ν])*(2.38^2)/6.
  SEIR_MCMC(25000, opt_cov, trace, priors, obs)
end

# Joint trace plots (last 100k iterations)
plotdf = DataFrame(iteration = rep(1:100000,6), value = [trace.α[end-99999:end], trace.β[end-99999:end], trace.ρ[end-99999:end], trace.γ[end-99999:end], trace.η[end-99999:end], trace.ν[end-99999:end]], parameter = [rep("α",100000),rep("β",100000),rep("η",100000),rep("ρ",100000),rep("γ",100000),rep("ν",100000)])
plot(plotdf, x="iteration", y="value", color="parameter", Geom.line)

# Histograms (last 100k iterations)
plot(x=trace.α[end-99999:end], Geom.histogram)
plot(x=trace.β[end-99999:end], Geom.histogram)
plot(x=trace.ρ[end-99999:end], Geom.histogram)
plot(x=trace.γ[end-99999:end], Geom.histogram)
plot(x=trace.η[end-99999:end], Geom.histogram)
