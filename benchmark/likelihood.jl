# Run with a Julia environment containing this checkout of Pathogen.
# Measurements exclude compilation and population/cache/workspace construction.
using Pathogen, Printf

constant_risk(p, pop, i) = p[1]
pair_risk(p, pop, i, j) = p[1]

function likelihood_fixture(n)
  pop = Population(DataFrame(x=zeros(n)))
  rf = RiskFunctions{SIR}(constant_risk, constant_risk, pair_risk, constant_risk, constant_risk)
  rp = RiskParameters{SIR}([0.01], [1.0], [1.0/n], [1.0], [0.5])
  k = n ÷ 2
  infection = [collect(0:k-1) .* 0.01; fill(NaN, n-k)]
  events = Events{SIR}(infection, infection .+ 0.05)
  return (rp, rf, events, pop, fill(State_S, n))
end

function measure_likelihood(n; reuse=true, cached=true)
  args = likelihood_fixture(n)
  cache = cached ? Pathogen.TransmissionRateCache(n) : nothing
  workspace = reuse ? Pathogen.LikelihoodWorkspace(SIR, n) : nothing
  evaluate() = Pathogen.loglikelihood(args...; transmission_rates_output=false,
    transmissions_output=false, transmission_rate_cache=cache, workspace=workspace)[1]
  evaluate()
  samples = [@timed evaluate() for _ in 1:5]
  @assert all(x -> isfinite(x[1]), samples)
  @printf("n=%d reuse=%s cached=%s ms=%.3f MiB=%.3f ll=%.12f\n", n, reuse, cached,
    1000minimum(x[2] for x in samples), minimum(x[3] for x in samples)/2.0^20, samples[1][1])
end

for n in (isempty(ARGS) ? [200, 400, 800, 1600] : parse.(Int, ARGS))
  for cached in (false, true)
    measure_likelihood(n; reuse=false, cached=cached)
    measure_likelihood(n; cached=cached)
  end
end
