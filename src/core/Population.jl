struct Population
  risks::AbstractDataFrame
  distances::Union{Nothing, AbstractArray}
  individuals::Int64

  Population(risks, distances) = new(risks, distances, size(risks, 1))
end

Population(risks) = Population(risks, nothing)

function Base.show(io::IO, x::Population)
  return print(io, "Population object (n=$(x.individuals))")
end
