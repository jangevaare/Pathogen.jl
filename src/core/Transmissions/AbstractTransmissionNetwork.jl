abstract type AbstractTransmissionNetwork end

function Base.sum(x::TN) where {TN <: AbstractTransmissionNetwork}
  return sum(x.external) + sum(x.internal)
end

function Base.show(io::IO, x::TN) where {TN <: AbstractTransmissionNetwork}
  return print(io, "$TN Σexternal = $(sum(x.external)) Σinternal = $(sum(x.internal))" )
end

function Base.copy(x::TN) where {TN <: AbstractTransmissionNetwork}
  return TN(copy(x.external),
            copy(x.internal))
end

function individuals(x::TN) where {TN <: AbstractTransmissionNetwork}
  return length(x.external)
end

function Base.sum(x::TN, i::Integer) where {TN <: AbstractTransmissionNetwork}
  if i < 1 || i > individuals(x)
    @error "Invalid individual identifier" i
  end
  return x.external[i] + sum(x.internal[:, i])
end

