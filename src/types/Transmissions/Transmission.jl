abstract type Transmission end

struct ExogenousTransmission
  individual::Int64
end

struct EndogenousTransmission
  individual::Int64
  source::Int64
end
