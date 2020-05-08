_subscript(i::Integer) = i<0 ? error("$i is negative") : join('₀'+d for d in reverse(digits(i)))

function _parameters(x::RiskParameters{M}) where {M <: DiseaseStateSequence}
  inds = _indices(x, cumulative=false)
  params1 = ["ϵ", "Ωs", "κ", "Ωt", "Ωl", "Ωr"]
  params2 = String[]
  for k = eachindex(inds)
    for l = 1:inds[k]
      push!(params2, params1[k] * _subscript(l))
    end
  end
  return params2
end

"""
Produces a data frame summarizing MCMC results.
"""
function Base.summary(mcmc::MCMC{M}; burnin::Int64=0, thin::Int64=1, bychain::Bool=true, credibleinterval::Float64=0.95) where M <: DiseaseStateSequence
  !(0.0 < credibleinterval < 1.0) && error("`credibleinterval` must be between 0 and 1")
  ci = (1-credibleinterval)/2
  params = _parameters(mcmc.markov_chains[1].risk_parameters[1])
  if bychain && length(mcmc.markov_chains) > 1
    df = DataFrame(parameter=String[],
                   chain=Int64[],
                   mean=Float64[],
                   var=Float64[],
                   CI=Tuple{Float64, Float64}[])
    for i = 1:length(mcmc.markov_chains)
      trace_array = convert(Array{Float64, 2}, mcmc.markov_chains[i].risk_parameters[1+burnin:thin:end])
      append!(df,
        DataFrame(parameter = params,
                  chain = fill(i, size(trace_array, 2)),
                  mean = mean(trace_array, dims=1)[:],
                  var = var(trace_array, dims=1)[:],
                  CI = [Tuple(quantile(trace_array[:,j], [ci, 1-ci])) for j = 1:size(trace_array, 2)]))
    end
  else
    trace_array = convert(Array{Float64, 2}, mcmc.markov_chains[1].risk_parameters[1+burnin:thin:end])
    for i = 2:length(mcmc.markov_chains)
      vcat!(trace_array, convert(Array{Float64, 2}, mcmc.markov_chains[i].risk_parameters[1+burnin:thin:end]))
    end
    df = DataFrame(parameter = params,
                   mean = mean(trace_array, dims=1)[:],
                   var = var(trace_array, dims=1)[:],
                   CI = [Tuple(quantile(trace_array[:,j], [ci, 1-ci])) for j = 1:size(trace_array, 2)])
  end
  return df
end