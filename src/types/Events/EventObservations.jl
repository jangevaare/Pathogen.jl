mutable struct EventObservations{T <: EpidemicModel}
  infection::Vector{Float64}
  removal::Vector{Float64}
  individuals::Int64

  function EventObservations{T}(i, r) where T <: Union{SEIR, SIR}
    if length(i) != length(r)
      @error "Length of infection and removal times must be equal"
    end
    return new{T}(i, r, length(i))
  end

  function EventObservations{T}(i) where T <: Union{SEI, SI}
    x = new{T}()
    x.infection = i
    x.individuals = length(i)
    return x
  end

  function EventObservations{T}(i::Array{Float64, 2}) where T <: Union{SEI, SI}
    if size(i, 2) != 1
      @error "Invalid Array dimensions for observations of a $T model"
    end
    x = new{T}()
    x.infection = i[:,1]
    x.individuals = size(i, 1)
    return x
  end

  function EventObservations{T}(ir::Array{Float64, 2}) where T <: Union{SEIR, SIR}
    if size(i, 2) != 2
      @error "Invalid Array dimensions for observations of a $T model"
    end
    x = new{T}()
    x.infection = ir[:, 1]
    x.removal = ir[:, 2]
    x.individuals = size(ir, 1)
    return x
  end
end

function Base.show(io::IO, x::EventObservations{T}) where T <: EpidemicModel
  return print(io, "$T model observations (n=$(x.individuals))")
end
