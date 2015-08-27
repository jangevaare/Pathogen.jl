"""
Example 2
SEIR simulation, visualization, and inference
Justin Angevaare
"""

using Pathogen, Gadfly, DataFrames, Distributions, ProgressMeter

cd("Desktop/example2")

# Simulate
init_seq = create_seq(400, 0.25, 0.25, 0.25, 0.25)
init_var = rand(Uniform(0,25), (2,400))

pop = create_population(init_seq, init_var)

powerlaw = create_powerlaw(3., 5., 0.001)
latency = create_constantrate(1/3.)
recovery = create_constantrate(1/5.)
substitution = jc69((0.01,))

ratearray = create_ratearray(pop, powerlaw, substitution)

n = 30000
progressbar = Progress(n, 5, "Simulating $n events...", 30)
@time while length(pop.timeline[1]) < n
  next!(progressbar)
  onestep!(ratearray, pop, powerlaw, latency, recovery, substitution)
end

# Inference
# α, β: powerlaw exposure kernel parameters
# η: external pressure rate
# ρ: infectivity rate (1/mean latent period)
# γ: recovery rate (1/mean infectious period)
# ν: detection rate (1/mean detection lag)

actual, obs = SEIR_surveilance(pop, 2.)

priors = SEIR_priors(Uniform(0,10), Uniform(0,10), Uniform(0,0.005), Uniform(0,1), Uniform(0,1), Uniform(1,3))

trace = SEIR_initialize(priors, obs)

SEIR_MCMC(100000, diagm([0.5, 0.5, 0.5, 0.5, 0.5, 0.5]), trace, priors, obs)

# Tune the transition kernel's covariance matrix over 300k iterations
n = 300
for i = 1:n

  # Progress bar
  if i == 1
    progressbar = Progress(n, 5, "Performing $n tuning MCMC stages...", 30)
  else
    next!(progressbar)
  end

  # Tune transition matrix
  opt_cov = cov([trace.α trace.β trace.ρ trace.γ trace.η trace.ν])*(2.38^2)/6.

  # Perform 1000 MCMC iterations
  SEIR_MCMC(1000, opt_cov, trace, priors, obs, false)
end

opt_cov = cov([trace.α trace.β trace.ρ trace.γ trace.η trace.ν])*(2.38^2)/6.
SEIR_MCMC(100000, opt_cov, trace, priors, obs)

# Simulation/Maximum posteriori visualization
images = 500
max_tracelp=findfirst(trace.logposterior.==maximum(trace.logposterior))

for time = 1:images
  states, routes = plotdata(pop, (time*maximum([maximum(trace.aug[max_tracelp]), pop.timeline[1][end]])/images))
  p1 = plot(layer(states, x="x", y="y", color="state", Geom.point),
            layer(routes, x="x", y="y", group="age", Geom.polygon),
            Theme(panel_opacity=1., panel_fill=color("white"), default_color=color("black"), background_color=color("white")))

  states, routes = plotdata(obs, trace, max_tracelp, (time*maximum([maximum(trace.aug[max_tracelp]), pop.timeline[1][end]])/images))
  p2 = plot(layer(states, x="x", y="y", color="state", Geom.point),
            layer(routes, x="x", y="y", group="age", Geom.polygon),
            Theme(panel_opacity=1., panel_fill=color("white"), default_color=color("black"), background_color=color("white")))

  filenumber = time/images
  filenumber = prod(split("$filenumber", ".", 2))
  filenumber *= prod(fill("0", 5-length(filenumber)))
  draw(PNG("SEIR_simulation_$filenumber.png", 15cm, 20cm), vstack(p1,p2))
end

# Assemble into animation
run(`convert -delay 10 -loop 0 -layers optimize SEIR_simulation_*.png SEIR_animation_combined.gif`)

# Remove frames
for time = 1:images
  filenumber = time/images
  filenumber = prod(split("$filenumber", ".", 2))
  filenumber *= prod(fill("0", 5-length(filenumber)))
  rm("SEIR_simulation_$filenumber.png")
end

# Inference visualization
# Joint trace plots (last 100k iterations)
plotdf = DataFrame(iteration = rep(1:100000,6), value = [trace.α[end-99999:end], trace.β[end-99999:end], trace.η[end-99999:end], trace.ρ[end-99999:end], trace.γ[end-99999:end], trace.ν[end-99999:end]], parameter = [rep("α",100000),rep("β",100000),rep("η",100000),rep("ρ",100000),rep("γ",100000),rep("ν",100000)])

draw(PNG("SEIR_traceplot.png", 20cm, 15cm),
     plot(plotdf,
          x="iteration",
          y="value",
          color="parameter",
          Geom.line,
          Theme(panel_opacity=1.,
                panel_fill=color("white"),
                background_color=color("white"))))

# logposterior plot (last 100k iterations)
draw(PNG("SEIR_logposterior.png", 20cm, 15cm),
     plot(x=1:100000,
          y=trace.logposterior[end-99999:end],
          Geom.line,
          Theme(panel_opacity=1.,
                panel_fill=color("white"),
                background_color=color("white"))))

# Posterior distribution histograms (last 100k iterations)
draw(PNG("SEIR_alpha_hist.png", 20cm, 15cm),
     plot(x=trace.α[end-99999:end],
          Geom.histogram,
          Theme(panel_opacity=1.,
                panel_fill=color("white"),
                background_color=color("white"))))

draw(PNG("SEIR_beta_hist.png", 20cm, 15cm),
     plot(x=trace.β[end-99999:end],
          Geom.histogram,
          Theme(panel_opacity=1.,
                panel_fill=color("white"),
                background_color=color("white"))))

draw(PNG("SEIR_eta_hist.png", 20cm, 15cm),
     plot(x=trace.η[end-99999:end],
          Geom.histogram,
          Theme(panel_opacity=1.,
                panel_fill=color("white"),
                background_color=color("white"))))

draw(PNG("SEIR_rho_hist.png", 20cm, 15cm),
     plot(x=trace.ρ[end-99999:end],
          Geom.histogram,
          Theme(panel_opacity=1.,
                panel_fill=color("white"),
                background_color=color("white"))))

draw(PNG("SEIR_gamma_hist.png", 20cm, 15cm),
     plot(x=trace.γ[end-99999:end],
          Geom.histogram,
          Theme(panel_opacity=1.,
                panel_fill=color("white"),
                background_color=color("white"))))

draw(PNG("SEIR_nu_hist.png", 20cm, 15cm),
     plot(x=trace.ν[end-99999:end],
          Geom.histogram,
          Theme(panel_opacity=1.,
                panel_fill=color("white"),
                background_color=color("white"))))

# Of those infected, what is the posterior probability of being exposed from external source (last 100k iterations)
network_posterior = mean(trace.network[end-100000:end])

y, x, z = findnz(network_posterior)
df = DataFrame(x=x, y=y, z=z)

draw(PNG("SEIR_exposure_network.png", 20cm, 20cm),
     plot(df,
          x="x",
          y="y",
          color="z",
          Scale.color_continuous(minvalue=0, maxvalue=1),
          Scale.x_continuous,
          Scale.y_continuous,
          Geom.rectbin,
          Stat.identity,
          Guide.xlabel(nothing),
          Guide.ylabel(nothing),
          Guide.colorkey("Posterior probability of exposure source"),
          Theme(panel_opacity=1.,
                panel_fill=color("white"),
                background_color=color("white"),
                key_position = :none)))
