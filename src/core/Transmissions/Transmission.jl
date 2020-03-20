abstract type Transmission end

struct ExogenousTransmission <: Transmission
  individual::Int64
end

struct EndogenousTransmission <: Transmission
  individual::Int64
  source::Int64
end

struct NoTransmission <: Transmission
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
