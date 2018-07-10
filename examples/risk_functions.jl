function _one(params::Vector{Float64}, pop::DataFrame, i::Int64)
  return 1.0
end

function _constant(params::Vector{Float64}, pop::DataFrame, i::Int64)
  return params[1]
end

function _coefficient(params::Vector{Float64}, pop::DataFrame, i::Int64)
  return params[1] * pop[:riskfactor1][i]
end

function _linear(params::Vector{Float64}, pop::DataFrame, i::Int64)
  return params[1] + params[2] * pop[:riskfactor1][i]
end

function _powerlaw(params::Vector{Float64}, pop::DataFrame, i::Int64, k::Int64)
  α = params[1]
  β = params[2]
  d = sqrt((pop[:x][i] - pop[:x][k])^2 + (pop[:y][i] - pop[:y][k])^2)
  return α * d^-β
end

function _powerlaw_w_intercept(params::Vector{Float64}, pop::DataFrame, i::Int64, k::Int64)
  α = params[1]
  β = params[2]
  γ = params[3]
  d = sqrt((pop[:x][i] - pop[:x][k])^2 + (pop[:y][i] - pop[:y][k])^2)
  return (α * d^-β) + γ
end

function _gaussian(params::Vector{Float64}, pop::DataFrame, i::Int64, k::Int64)
  σ = params[1]
  d = sqrt((pop[:x][i] - pop[:x][k])^2 + (pop[:y][i] - pop[:y][k])^2)/σ
  return pdf(Normal(), d)
end

function _gaussian_w_intercept(params::Vector{Float64}, pop::DataFrame, i::Int64, k::Int64)
  σ = params[1]
  γ = params[2]
  d = sqrt((pop[:x][i] - pop[:x][k])^2 + (pop[:y][i] - pop[:y][k])^2)/σ
  return pdf(Normal(), d) + γ
end
