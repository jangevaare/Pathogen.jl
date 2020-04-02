function generate(::Type{SM},
                  priors::Vector{UnivariateDistribution}) where {
                  SM <: NucleicAcidSubstitutionModel}
  return SM([rand(x) for x in priors])
end


function generate(::Type{SM},
                  last_sm::SM,
                  Σ::Array{Float64, 2}) where {
                  SM <: NucleicAcidSubstitutionModel}
  if size(Σ) == (0,0)
    return last_sm
  else
    return SM(rand(MvNormal([getproperty(last_sm, θ) for θ in [propertynames(SM)...]], Σ)), safe=false)
  end
end