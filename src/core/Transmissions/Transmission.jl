abstract type AbstractTransmission end

struct ExogenousTransmission <: AbstractTransmission
  individual::Int64
end

struct EndogenousTransmission <: AbstractTransmission
  individual::Int64
  source::Int64
end

struct NoTransmission <: AbstractTransmission
end

function Base.show(io::IO, x::ExogenousTransmission)
  return println(io, "Exogenous infection transmission to individual $(x.individual)")
end

function Base.show(io::IO, x::EndogenousTransmission)
  return println(io, "Endogenous infection transmission to individual $(x.individual) by individual $(x.source)")
end

function Base.show(io::IO, x::NoTransmission)
  return println(io, "No infection transmission")
end
