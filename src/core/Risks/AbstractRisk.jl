abstract type AbstractRisk{S <: DiseaseStateSequence} end

function Base.copy(x::T) where {S <: SEIR, T <: AbstractRisk{S}}
  return T{S}(copy(x.sparks),
              copy(x.susceptibility),
              copy(x.infectivity),
              copy(x.transmissibility),
              copy(x.latency),
              copy(x.removal))
end

function Base.copy(x::T) where {S <: SEI, T <: AbstractRisk{S}}
  return T{S}(copy(x.sparks),
              copy(x.susceptibility),
              copy(x.infectivity),
              copy(x.transmissibility),
              copy(x.latency))
end

function Base.copy(x::T) where {S <: SIR, T <: AbstractRisk{S}}
  return T{S}(copy(x.sparks),
              copy(x.susceptibility),
              copy(x.infectivity),
              copy(x.transmissibility),
              copy(x.removal))
end

function Base.copy(x::T) where {S <: SI, T <: AbstractRisk{S}}
  return T{S}(copy(x.sparks),
              copy(x.susceptibility),
              copy(x.infectivity),
              copy(x.transmissibility))
end

function _indices(x::T; zeros::Bool=true) where {S <: DiseaseStateSequence, T <: AbstractRisk{S}}
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


function Base.length(x::T) where {S <: DiseaseStateSequence, T <: AbstractRisk{S}}
  return _indices(x)[end]
end


function Base.getindex(x::T, i::Int64) where {S <: DiseaseStateSequence, T <: AbstractRisk{S}}
  indices = _indices(x, zeros = true)
  riskfunc = findfirst(i .<= indices)
  return getfield(x, riskfunc)[end - (indices[riskfunc] - i)]
end
