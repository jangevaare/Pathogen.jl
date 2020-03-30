primitive type DiseaseState 8 end

DiseaseState(x::UInt8) = reinterpret(DiseaseState, x)
UInt8(x::DiseaseState) = reinterpret(UInt8, x)

# `DieaseState` const representations
const State_S = DiseaseState(0b0001)
const State_E = DiseaseState(0b0010)
const State_I = DiseaseState(0b0100)
const State_R = DiseaseState(0b1000)

const DiseaseStates = Vector{DiseaseState}

const _DiseaseStateVector = SVector{4}(State_S, State_E, State_I, State_R)
const _DiseaseStateCharVector = SVector{4}('S', 'E', 'I', 'R')

# `DiseaseState` conversion to/from `Char`
function Base.convert(::Type{Char}, x::DiseaseState)
  return _DiseaseStateCharVector[_DiseaseStateVector .== Ref(x)][1]
end

function Base.show(io::IO, x::DiseaseState)
  return print(io, convert(Char, x))
end

Base.convert(::Type{DiseaseStates}, x::Type{SEIR}) = SVector{4}(State_S, State_E, State_I, State_R)

function Base.convert(::Type{DiseaseStates}, x::Type{SEI})
  return SVector{3}(State_S, State_E, State_I)
end

function Base.convert(::Type{DiseaseStates}, x::Type{SIR})
  return SVector{3}(State_S, State_I, State_R)
end

function Base.convert(::Type{DiseaseStates}, x::Type{SI})
  return SVector{2}(State_S, State_I)
end