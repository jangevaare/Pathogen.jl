function mean(x::Vector{Events{T}}) where T <: EpidemicModel
  return Events{T}(mean([convert(Array{Float64, 2}, i) for i in x]))
end
