function generate(::Type{SM},
                  priors::Vector{UnivariateDistribution}) where {
                  SM <: NucleicAcidSubstitutionModel}
  return SM([rand(x) for x in priors])
end


function generate(::Type{SM},
                  last_sm::SM,
                  Σ::Array{FLoat64, 2}, where {
                  SM <: NucleicAcidSubstitutionModel}
 return SM(rand(MvNormal([getproperty(last_sm, θ) for θ in [propertynames(x1)...]], Σ)), safe=false)
end