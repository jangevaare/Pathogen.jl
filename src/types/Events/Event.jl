struct Event{T <: EpidemicModel}
  time::Float64
  individual::Int64
  prior::DiseaseState
  present::DiseaseState
  # TODO Add `EpidemicModel` specific `DiseaseState` validation
end

# abstract type Event{T <: EpidemicModel} end
#
# struct Exposure{T <: Union{SEIR, SEI}} <: Event{T}
#   individual::Int64
# end
#
# struct Infection{T <: Union{SEIR, SEI, SIR, SI}} <: Event{T}
#   individual::Int64
# end
#
# struct Removal{T <: Union{SEIR, SIR}} <: Event{T}
#   individual::Int64
# end
