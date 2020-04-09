_subscript(i::Integer) = i<0 ? error("$i is negative") : join('₀'+d for d in reverse(digits(i)))

"""
Produces a data frame summarizing MCMC results.
"""
function summary(mcmc::MCMC{M}; burnin::Int64=0, thin::Int64=1, bychain::Bool=true) where M <: EpidemicModel
  if bychain && length(mcmc.markov_chains) > 1
    df = DataFrame(parameter=String[],
                   chain=Int64[],
                   mean=Float64[],
                   var=Float64[],
                   credible99=Tuple{Float64, Float64}[],
                   credible95=Tuple{Float64, Float64}[],
                   credible90=Tuple{Float64, Float64}[])
    for i = 1:length(mcmc.markov_chains)
      trace_array = convert(Array{Float64, 2}, mcmc.markov_chains[i].risk_parameters[1+burnin:thin:end])
      inds = _indices(mcmc.markov_chains[i].risk_parameters[1], cumulative=false)
      params1 = ["ϵ", "Ωs", "κ", "Ωt", "Ωl", "Ωr"]
      params2 = String[]
      for k = eachindex(inds)
        for l = 1:inds[k]
          push!(params2, params1[k] * _subscript(l))
        end
      end
      append!(df,
        DataFrame(parameter = params2,
                  chain = fill(i, size(trace_array, 2)),
                  mean = mean(trace_array, dims=1)[:],
                  var = var(trace_array, dims=1)[:],
                  credible99 = [Tuple(quantile(trace_array[:,j], [0.005, 0.995])) for j = 1:size(trace_array, 2)],
                  credible95 = [Tuple(quantile(trace_array[:,j], [0.025, 0.975])) for j = 1:size(trace_array, 2)],
                  credible90 = [Tuple(quantile(trace_array[:,j], [0.05, 0.95])) for j = 1:size(trace_array, 2)]))
    end
  else
    trace_array = convert(Array{Float64, 2}, mcmc.markov_chains[1].risk_parameters[1+burnin:thin:end])
    for i = 2:length(mcmc.markov_chains)
      vcat!(trace_array, convert(Array{Float64, 2}, mcmc.markov_chains[i].risk_parameters[1+burnin:thin:end]))
    end
    inds = _indices(mcmc.markov_chains[1].risk_parameters[1], cumulative=false)
    params1 = ["ϵ", "Ωs", "κ", "Ωt", "Ωl", "Ωr"]
    params2 = String[]
    for k = eachindex(inds)
      for l = 1:inds[k]
        push!(params2, params1[k] * _subscript(l))
      end
    end
    df = DataFrame(parameter = params2,
                   mean = mean(trace_array, dims=1)[:],
                   var = var(trace_array, dims=1)[:],
                   credible99 = [Tuple(quantile(trace_array[:,j], [0.005, 0.995])) for j = 1:size(trace_array, 2)],
                   credible95 = [Tuple(quantile(trace_array[:,j], [0.025, 0.975])) for j = 1:size(trace_array, 2)],
                   credible90 = [Tuple(quantile(trace_array[:,j], [0.050, 0.950])) for j = 1:size(trace_array, 2)])
  end
  return df
end