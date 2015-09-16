"""
Example 1
SEIR simulation, visualization, and inference
Justin Angevaare
"""
using Pathogen, Gadfly, DataFrames, Distributions, ProgressMeter

cd("Desktop/example1")

# Simulate
init_seq = create_seq(100, 0.25, 0.25, 0.25, 0.25)
init_var = rand(Uniform(0,10), (2,40))

pop = create_population(init_seq, init_var)

powerlaw = create_powerlaw(4., 5., 0.001)
latency = create_constantrate(1/3.)
recovery = create_constantrate(1/5.)
substitution = jc69([0.05])

ratearray = create_ratearray(pop, powerlaw, substitution)

@time while length(pop.timeline[1]) < 5000.
  onestep!(ratearray, pop, powerlaw, latency, recovery, substitution)
end

# Inference
# α, β: powerlaw exposure kernel parameters
# η: external pressure rate
# ρ: infectivity rate (1/mean latent period)
# γ: recovery rate (1/mean infectious period)
# ν: detection rate (1/mean detection lag)

actual, obs = surveil(pop, 2.)

ilm_priors = SEIR_priors(Uniform(2,8), Uniform(2,8), Uniform(0,0.005), Uniform(0,1), Uniform(0,1))
detection_priors = Lag_priors(Uniform(1,3))
mutation_priors = JC69_priors(Uniform(0,0.1))

# From some previous simulation runs...
opt_cov = [0.23086898731119626 -0.06339624871696066 -0.0031800042117832114 0.01611130363856698 7.502308677590557e-5 -0.04172552030369371 9.656673691052941e-6
 -0.06339624871696066 0.7402072386878623 -0.06563597504892492 -0.020217534870576297 -0.0002459116200602692 0.1563007257329955 -5.731936222295633e-6
 -0.0031800042117832114 -0.06563597504892492 0.036204426712783414 0.00014101116046161534 3.31375867657149e-5 -0.019469166893799742 2.2529637404798412e-6
 0.01611130363856698 -0.020217534870576297 0.00014101116046161534 0.009196894826401476 2.7716042390333208e-5 -0.00574364739285761 5.448990926026132e-7
 7.502308677590557e-5 -0.0002459116200602692 3.31375867657149e-5 2.7716042390333208e-5 1.7177313965323465e-6 -2.4717221007601892e-5 7.801292516978605e-9
 -0.04172552030369371 0.1563007257329955 -0.019469166893799742 -0.00574364739285761 -2.4717221007601892e-5 0.15638126601298286 -5.625526286964686e-6
 9.656673691052941e-6 -5.731936222295633e-6 2.2529637404798412e-6 5.448990926026132e-7 7.801292516978605e-9 -5.625526286964686e-6 1.6749597385278377e-8]

ilm_trace, detection_trace, mutation_trace = MCMC(100000, opt_cov, ilm_priors, detection_priors, mutation_priors, obs, false, true, true)

# Tune the transition kernel's covariance matrix
n = 100
for i = 1:n

  # Progress bar
  if i == 1
    progressbar = Progress(n, 5, "Performing $n tuning MCMC stages...", 30)
  else
    next!(progressbar)
  end

  # Tune transition matrix
  opt_cov = cov([ilm_trace.α ilm_trace.β ilm_trace.ρ ilm_trace.γ ilm_trace.η detection_trace.ν mutation_trace.λ])*(2.38^2)/7.

  # Perform 1000 MCMC iterations
  MCMC(1000, opt_cov, ilm_trace, detection_trace, mutation_trace, ilm_priors, detection_priors, mutation_priors, obs, false, false, false)
end

opt_cov = cov([ilm_trace.α ilm_trace.β ilm_trace.ρ ilm_trace.γ ilm_trace.η detection_trace.ν mutation_trace.λ])*(2.38^2)/7.

MCMC(100000, opt_cov, ilm_trace, detection_trace, mutation_trace, ilm_priors, detection_priors, mutation_priors, obs, false, true, true)

# Simulation/Maximum posteriori visualization
images = 500
max_tracelp=findfirst(ilm_trace.logposterior.==maximum(ilm_trace.logposterior))

for time = 1:images
  states, routes = plotdata(pop, (time*maximum([maximum(ilm_trace.aug[max_tracelp]), pop.timeline[1][end]])/images))
  p1 = plot(layer(states, x="x", y="y", color="state", Geom.point),
            layer(routes, x="x", y="y", group="age", Geom.polygon),
            Theme(panel_opacity=1., panel_fill=color("white"), default_color=color("black"), background_color=color("white")))

  states, routes = plotdata(obs, ilm_trace, max_tracelp, (time*maximum([maximum(ilm_trace.aug[max_tracelp]), pop.timeline[1][end]])/images))
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
plotdf = DataFrame(iteration = rep(1:100000,7),
                   value = [ilm_trace.α[end-99999:end],
                            ilm_trace.β[end-99999:end],
                            ilm_trace.η[end-99999:end],
                            ilm_trace.ρ[end-99999:end],
                            ilm_trace.γ[end-99999:end],
                            detection_trace.ν[end-99999:end],
                            mutation_trace.λ[end-99999:end]],
                   parameter = [rep("α",100000),
                                rep("β",100000),
                                rep("η",100000),
                                rep("ρ",100000),
                                rep("γ",100000),
                                rep("ν",100000),
                                rep("λ",100000)])

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
          y=ilm_trace.logposterior[end-99999:end],
          Geom.line,
          Theme(panel_opacity=1.,
                panel_fill=color("white"),
                background_color=color("white"))))

# Posterior distribution histograms (last 100k iterations)
draw(PNG("SEIR_alpha_hist.png", 20cm, 15cm),
     plot(x=ilm_trace.α[end-99999:end],
          Geom.histogram,
          Theme(panel_opacity=1.,
                panel_fill=color("white"),
                background_color=color("white"))))

draw(PNG("SEIR_beta_hist.png", 20cm, 15cm),
     plot(x=ilm_trace.β[end-99999:end],
          Geom.histogram,
          Theme(panel_opacity=1.,
                panel_fill=color("white"),
                background_color=color("white"))))

draw(PNG("SEIR_eta_hist.png", 20cm, 15cm),
     plot(x=ilm_trace.η[end-99999:end],
          Geom.histogram,
          Theme(panel_opacity=1.,
                panel_fill=color("white"),
                background_color=color("white"))))

draw(PNG("SEIR_rho_hist.png", 20cm, 15cm),
     plot(x=ilm_trace.ρ[end-99999:end],
          Geom.histogram,
          Theme(panel_opacity=1.,
                panel_fill=color("white"),
                background_color=color("white"))))

draw(PNG("SEIR_gamma_hist.png", 20cm, 15cm),
     plot(x=ilm_trace.γ[end-99999:end],
          Geom.histogram,
          Theme(panel_opacity=1.,
                panel_fill=color("white"),
                background_color=color("white"))))

draw(PNG("SEIR_nu_hist.png", 20cm, 15cm),
     plot(x=detection_trace.ν[end-99999:end],
          Geom.histogram,
          Theme(panel_opacity=1.,
                panel_fill=color("white"),
                background_color=color("white"))))

draw(PNG("SEIR_lambda_hist.png", 20cm, 15cm),
     plot(x=mutation_trace.λ[end-99999:end],
          Geom.histogram,
          Theme(panel_opacity=1.,
                panel_fill=color("white"),
                background_color=color("white"))))

# Of those infected, what is the posterior probability of being exposed from external source (last 100k iterations)
network_posterior = mean(ilm_trace.network[end-100000:end])

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
