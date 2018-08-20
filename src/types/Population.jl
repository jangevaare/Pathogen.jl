mutable struct Population
  risks::DataFrame
  distances::Array{Any, 2}
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

  function Population(d::Array{Any, 2})
    if length(unique(size(d))) != 1
      @error "Distance matrices must be square"
    end
    x = new()
    x.distances = d
    x.individuals = size(d, 1)
    return x
  end

  function Population(risks::DataFrame, d::Array{Float64, 2})
    if length(unique([size(d, 1); size(d, 2); size(risks, 1)])) !== 1
      @error "Mismatch between sizes of distance matrix and risk dataframe"
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
