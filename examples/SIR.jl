using Distances,
      DataFrames,
      Distributions,
      Pathogen

n = 100
risks = DataFrame(x = rand(Uniform(0, 20), n),
                  y = rand(Uniform(0, 40), n),
                  riskfactor1 = rand(Gamma(), n))

# Will precalculate distances
dists = [euclidean([risks[i, :x];
                    risks[i, :y]],
                   [risks[j, :x];
                    risks[j, :y]]) for i = 1:n, j = 1:n]

pop = Population(risks, dists)

function _constant(params::Vector{Float64}, pop::Population, i::Int64)
  return params[1]
end

function _one(params::Vector{Float64}, pop::Population, i::Int64)
  return 1.0
end

function _linear(params::Vector{Float64}, pop::Population, i::Int64)
  return params[1] * pop.risks[i, :riskfactor1]
end

function _powerlaw(params::Vector{Float64}, pop::Population, i::Int64, k::Int64)
  α = params[1]
  β = params[2]
  d = pop.distances[k, i]
  return α * (d^(-β))
end

rf = RiskFunctions{SIR}(_constant, # sparks function
                        _one, # susceptibility function
                        _powerlaw, # infectivity function
                        _one, # transmissability function
                        _linear) # removal function

rparams = RiskParameters{SIR}([0.0001], # sparks function parameter(s)
Float64[], # susceptibility function parameter(s)
[3.0, 4.0], # infectivity function parameter(s)
Float64[], # transmissibility function parameter(s)
[0.05]) # removal function parameter(s)

starting_states = append!([State_I], fill(State_S, n-1)) # Set first individual as infectious, others as susceptible to start

sim = Simulation(pop, starting_states, rf, rparams)

simulate!(sim, tmax=100.0)

using Plots, Plots.PlotMeasures
#pgfplots()
gr(dpi=400)
#unicodeplots()


p0 = plot(sim.events, 0.0, 100.0, linewidth=2.0, legendfont=font(6), xaxis=font(10), bottom_margin=30px)
p1=plot(sim.transmission_network, sim.population, sim.events, 0.0, linealpha=1.0, linewidth=1.0, linecolour=:black, axis=nothing, foreground_color_subplot=:white, legend=:none, markersize=2.75, markerstrokewidth=1.0, markerstrokealpha=1.0, markerstrokecolour=:black, title="Time = 0", titlefontcolor=:black, titlefontsize = 8)
p2 =plot(sim.transmission_network, sim.population, sim.events, 10.0, linealpha=1.0, linewidth=1.0, linecolour=:black, axis=nothing, foreground_color_subplot=:white, legend=:none, markersize=2.75, markerstrokewidth=1.0, markerstrokealpha=1.0, markerstrokecolour=:black, title="Time = 10", titlefontcolor=:black, titlefontsize = 8)
p3 =plot(sim.transmission_network, sim.population, sim.events, 30.0, linealpha=1.0, linewidth=1.0, linecolour=:black, axis=nothing, foreground_color_subplot=:white, legend=:none, markersize=2.75, markerstrokewidth=1.0, markerstrokealpha=1.0, markerstrokecolour=:black, title="Time = 30", titlefontcolor=:black, titlefontsize = 8)
p4 =plot(sim.transmission_network, sim.population, sim.events, 60.0, linealpha=1.0, linewidth=1.0, linecolour=:black, axis=nothing, foreground_color_subplot=:white, legend=:none, markersize=2.75, markerstrokewidth=1.0, markerstrokealpha=1.0, markerstrokecolour=:black, title="Time = 60", titlefontcolor=:black, titlefontsize = 8)
p5 =plot(sim.transmission_network, sim.population, sim.events, 100.0, linealpha=1.0, linewidth=1.0, linecolour=:black, axis=nothing, foreground_color_subplot=:white, legend=:none, markersize=2.75, markerstrokewidth=1.0, markerstrokealpha=1.0, markerstrokecolour=:black, title="Time = 100", titlefontcolor=:black, titlefontsize = 8)
l1 = @layout [a;
              b c d e f]
r1 = plot(p0, p1, p2, p3, p4, p5, layout=l1)
png(r1, "epiplot.png")

# Generate observations with Uniform(0, 2) observation delay for infection and removal
obs = observe(sim, Uniform(0.0, 2.0), Uniform(0.0, 2.0), force=true)

# Optimistically assume we know the functional form of epidemic (i.e. use same risk functions used for simulation purposes)
# Specify some priors for the risk parameters of our various risk functions
# Set some extents for event data augmentation

rpriors = RiskPriors{SIR}([Exponential(0.0001)],
                          UnivariateDistribution[],
                          [Uniform(0.0, 2.0); Uniform(1.0, 8.0)],
                          UnivariateDistribution[],
                          [Uniform(0.0, 1.0)])

ee = EventExtents{SIR}(5.0, 5.0)

# Initialize MCMC
mcmc = MCMC(obs, ee, pop, rf, rpriors)
start!(mcmc, attempts=25000) # 1 chain, with 10k initialization attempts each

# Run MCMC
start_time = time()
iterate!(mcmc, 25000, 1.0, condition_on_network=true, event_batches=10)
mcmc_time = time() - start_time

p0 = plot(1:20:25001, convert(Array{Float64,2}, mcmc.markov_chains[1].risk_parameters[1:20:25001]), yscale=:log10, legend=:none, title="TN-ILM parameters", xlab="Iteration", ylab="Value", xguidefontsize=8, yguidefontsize=8, xtickfontsize=7, ytickfontsize=7, titlefontsize=11, bottom_margin=30px)

p1 = plot(Pathogen._epidemic_curve(mcmc.markov_chains[1].events[5000], State_S, 0.0, 100.0),
          alpha=0.0125, legend=:none, colour=1,
          dpi=200, linewidth=1.0, seriestype=:steppost,
          ylim=(0,n), ylab="N",
          xlim=(0,100), xlab="Time",
          title="S", xguidefontsize=8, yguidefontsize=8, xtickfontsize=7, ytickfontsize=7, titlefontsize=11)
for i=5020:20:25000
plot!(p1, Pathogen._epidemic_curve(mcmc.markov_chains[1].events[i], State_S, 0.0, 100.0), alpha=0.0125, colour=1, linewidth=2.0, seriestype=:steppost)
end
plot!(p1, Pathogen._epidemic_curve(sim.events, State_S, 0.0, 100.0), colour=:black, linewidth=2.0, seriestype=:steppost)

p2= plot(Pathogen._epidemic_curve(mcmc.markov_chains[1].events[5000], State_I, 0.0, 100.0),
         alpha=0.0125, legend=:none, colour=1,
         dpi=400, linewidth=1.0, seriestype=:steppost,
         ylim=(0,n), ylab="N",
         xlim=(0,100), xlab="Time",
         title="I", xguidefontsize=8, yguidefontsize=8, xtickfontsize=7, ytickfontsize=7, titlefontsize=11)

for i=5020:20:25000
plot!(p2, Pathogen._epidemic_curve(mcmc.markov_chains[1].events[i], State_I, 0.0, 100.0), alpha=0.0125, colour=1, linewidth=2.0, seriestype=:steppost)
end
plot!(p2, Pathogen._epidemic_curve(sim.events, State_I, 0.0, 100.0), colour=:black, linewidth=2.0, seriestype=:steppost)

p3 = plot(Pathogen._epidemic_curve(mcmc.markov_chains[1].events[5000], State_R, 0.0, 100.0),
          alpha=0.0125, legend=:none, colour=1,
          dpi=400, linewidth=1.0, seriestype=:steppost,
          ylim=(0,n), ylab="N",
          xlim=(0,100), xlab="Time",
          title="R", xguidefontsize=8, yguidefontsize=8, xtickfontsize=7, ytickfontsize=7, titlefontsize=11)
for i=5020:20:25000
  plot!(p3, Pathogen._epidemic_curve(mcmc.markov_chains[1].events[i], State_R, 0.0, 100.0), alpha=0.0125, colour=1, linewidth=2.0, seriestype=:steppost)
end
plot!(p3, Pathogen._epidemic_curve(sim.events, State_R, 0.0, 100.0), colour=:black, linewidth=2.0, seriestype=:steppost)

l = @layout [a; [b c d]]
plot(p0, p1, p2, p3, layout=l)

png("posterior.png")

p4=plot(sim.transmission_network, sim.population, sim.events, 100.0, colour=:black, legend=:none, linealpha=1.0, linewidth=2.0, linecolour=:black, axis=nothing, foreground_color_subplot=:white)
p5=plot(mcmc.markov_chains[1].transmission_network[5020], sim.population, sim.events, 100.0, colour=:black, legend=:none, linealpha=0.005, axis=nothing, foreground_color_subplot=:white)
for i=5000:50:25000
  plot!(mcmc.markov_chains[1].transmission_network[i], sim.population, sim.events, 100.0, colour=:black, legend=:none, linealpha=0.005)
end
l = @layout [a b]
plot(p4, p5, layout=l)
png("posterior_tn_sbs.png")
