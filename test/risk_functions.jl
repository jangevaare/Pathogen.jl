function _one(params::Vector{Float64}, pop::Pathogen.Population, i::Int64)
  return 1.0
end

function _zero(params::Vector{Float64}, pop::Pathogen.Population, i::Int64)
  return 0.0
end

function _constant(params::Vector{Float64}, pop::Pathogen.Population, i::Int64)
  return params[1]
end

function _coefficient(params::Vector{Float64}, pop::Pathogen.Population, i::Int64)
  return params[1] * pop.risks[i, :riskfactor1]
end

function _linear(params::Vector{Float64}, pop::Pathogen.Population, i::Int64)
  return params[1] + params[2] * pop.risks[i, :riskfactor1]
end

function _powerlaw(params::Vector{Float64}, pop::Pathogen.Population, i::Int64, k::Int64)
  α = params[1]
  β = params[2]
  d = pop.distances[k, i]
  return α * (d^(-β))
end

function _powerlaw_w_intercept(params::Vector{Float64}, pop::Pathogen.Population, i::Int64, k::Int64)
  α = params[1]
  β = params[2]
  γ = params[3]
  d = pop.distances[k, i]
  return (α * (d^(-β))) + γ
end
