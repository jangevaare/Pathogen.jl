mutable struct EventObservations{T <: EpidemicModel}
  infection::Vector{Float64}
  removal::Vector{Float64}
  individuals::Int64

  function EventObservations{T}(i, r) where T <: Union{SEIR, SIR}
    if length(i) != length(r)
      error("Length of infection and removal times must be equal")
    end
    return new{T}(i, r, length(i))
  end

  function EventObservations{T}(i) where T <: Union{SEI, SI}
    x = new{T}()
    x.infection = i
    x.individuals = length(i)
    return x
  end
end
