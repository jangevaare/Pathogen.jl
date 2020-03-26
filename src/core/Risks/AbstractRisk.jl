abstract type AbstractRisk{S <: DiseaseStateSequence} end

function Base.copy(x::AbstractRisk{S}) where {S <: SEIR}
  return AbstractRisk{S}(copy(x.sparks),
                         copy(x.susceptibility),
                         copy(x.infectivity),
                         copy(x.transmissibility),
                         copy(x.latency),
                         copy(x.removal))
end

function Base.copy(x::AbstractRisk{S}) where {S <: SEI}
  return AbstractRisk{S}(copy(x.sparks),
                         copy(x.susceptibility),
                         copy(x.infectivity),
                         copy(x.transmissibility),
                         copy(x.latency))
end

function Base.copy(x::AbstractRisk{S}) where {S <: SIR}
  return AbstractRisk{S}(copy(x.sparks),
                         copy(x.susceptibility),
                         copy(x.infectivity),
                         copy(x.transmissibility),
                         copy(x.removal))
end

function Base.copy(x::AbstractRisk{S}) where {S <: SI}
  return AbstractRisk{S}(copy(x.sparks),
                         copy(x.susceptibility),
                         copy(x.infectivity),
                         copy(x.transmissibility))
end

function _indices(x::AbstractRisk{S}; zeros::Bool=true) where {S <: DiseaseStateSequence}
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
  return cumsum(indices)
end


function Base.length(x::AbstractRisk{S}) where {S <: DiseaseStateSequence}
  return _indices(x)[end]
end


function Base.getindex(x::AbstractRisk{S}, i::Int64) where {S <: DiseaseStateSequence}
  indices = _indices(x, zeros = true)
  riskfunc = findfirst(i .<= indices)
  return getfield(x, riskfunc)[end - (indices[riskfunc] - i)]
end
