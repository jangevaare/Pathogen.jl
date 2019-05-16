mutable struct Population
  risks::DataFrame
  distances::AbstractArray{T,2} where T <: Any
  individuals::Int64

  function Population(n::Int64)
    x = new()
    x.individuals = n
    return x
  end

  function Population(risks::DataFrame)
    x = new()
    x.risks = risks
    x.individuals = size(risks, 1)
    return x
  end

  function Population(d::AbstractArray{T,2}) where T <: Any
    if size(d, 1) !== size(d, 2)
      @error "Distance matrix must be square"
    end
    x = new()
    x.distances = d
    x.individuals = size(d, 1)
    return x
  end

  function Population(risks::DataFrame, d::AbstractArray{T,2}) where T <: Any
    if size(d, 1) !== size(d, 2)
      @error "Distance matrix must be square"
    elseif size(d, 1) !== size(risks, 1)
      @error "Mismatch between size of distance matrix and risk dataframe"
    end
    x = new()
    x.risks = risks
    x.distances = d
    x.individuals = size(d, 1)
    return x
  end
end

function Base.show(io::IO, x::Population)
  return print(io, "Population object (n=$(x.individuals)")
end
