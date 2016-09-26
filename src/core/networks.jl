"""
An infectious disease transmission network
"""
type Network
  external::Vector{Bool}
  internal::Array{Bool, 2}
  function Network(population::DataFrame)
    individuals = size(population, 1)
    return new(fill(false, individuals),
               fill(false, (individuals, individuals)))
  end

  function Network(external, internal)
    if !(length(external) == size(internal, 1) == size(internal, 2))
      error("Dimension mismatch in exposure network object")
    end
    return new(external, internal)
  end
end
