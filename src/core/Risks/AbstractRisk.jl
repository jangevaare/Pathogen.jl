abstract type AbstractRisk{S <: DiseaseStateSequence} end

function _indices(x::AbstractRisk{S};
                  zeros::Bool=true,
                  cumulative::Bool=true) where {
                  S <: DiseaseStateSequence}
  indices = [length(x.sparks)
             length(x.susceptibility)
             length(x.infectivity)
             length(x.transmissibility)]
  if S in [SEIR; SEI]
    push!(indices, length(x.latency))
  elseif zeros
    push!(indices, 0)
  end
  if S in [SEIR; SIR]
    push!(indices, length(x.removal))
  elseif zeros
    push!(indices, 0)
  end
  if cumulative
    return cumsum(indices)
  else
    return indices
  end
end

function Base.length(x::T) where {S <: DiseaseStateSequence, T <: AbstractRisk{S}}
  return _indices(x)[end]
end

function Base.getindex(x::T, i::Int64) where {S <: DiseaseStateSequence, T <: AbstractRisk{S}}
  indices = _indices(x, zeros = true)
  riskfunc = findfirst(i .<= indices)
  return getfield(x, riskfunc)[end - (indices[riskfunc] - i)]
end
