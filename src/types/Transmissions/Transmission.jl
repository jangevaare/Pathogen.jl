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
