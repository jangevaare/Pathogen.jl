struct Population{R <: AbstractDataFrame, D <: Union{Nothing, AbstractArray}}
  risks::R
  distances::D
end

Population(risks) = Population(risks, nothing)

function individuals(x::Population)
  return size(x.risks, 1)
end

function Base.show(io::IO, x::Population)
  return print(io, "Population object (n=$(individuals(x)))")
end