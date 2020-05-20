@recipe function f(x::Vector{RiskParameters{T}}) where T <: DiseaseStateSequence
  xguide --> "Iteration"
  yguide --> "Value"
  legend --> :none
  # label --> permutedims(_parameters(x[1]))
  convert(Array{Float64, 2}, x)
end

@recipe function f(iter,
                   x::Vector{RiskParameters{T}}) where T <: DiseaseStateSequence
  xguide --> "Iteration"
  yguide --> "Value"
  legend --> :none
  # label --> permutedims(_parameters(x[1]))
  iter, convert(Array{Float64, 2}, x[iter])
end

@recipe function f(x::Vector{T}) where {T <: NASM}
  xguide --> "Iteration"
  yguide --> "Value"
  legend --> :none
  # label --> permutedims(_parameters(x[1]))
  convert(Array{Float64, 2}, x)
end

@recipe function f(iter, x::Vector{T}) where {T <: NASM}
  xguide --> "Iteration"
  yguide --> "Value"
  legend --> :none
  # label --> permutedims(_parameters(x[1]))
  iter, convert(Array{Float64, 2}, x[iter])
end
