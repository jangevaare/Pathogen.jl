primitive type DiseaseState 8 end

DiseaseState(x::UInt8) = reinterpret(DiseaseState, x)
UInt8(x::DiseaseState) = reinterpret(UInt8, x)

# `DieaseState` const representations
const State_S = DiseaseState(0b0001)
const State_E = DiseaseState(0b0010)
const State_I = DiseaseState(0b0100)
const State_R = DiseaseState(0b1000)

const _DiseaseStateVector = [State_S; State_E; State_I; State_R]
const _DiseaseStateCharVector = ['S'; 'E'; 'I'; 'R']

# `DiseaseState` conversion to/from `Char`
function Base.convert(::Type{Char}, x::DiseaseState)
  return _DiseaseStateCharVector[_DiseaseStateVector .== x][1]
end

function Base.convert(::Type{DiseaseState}, x::Char)
  return _DiseaseStateVector[_DiseaseStateCharVector .== x][1]
end

# `DiseaseState` conversion to/from `BitArray`
function Base.convert(::Type{BitArray{1}}, x::DiseaseState)
  return BitArray(_DiseaseStateVector .== x)
end

function Base.convert(::Type{DiseaseState}, x::BitArray{1})
  return _DiseaseStateVector[x]
end

# `DiseaseState` conversion to/from `Int64`
function Base.convert(::Type{Int64}, x::DiseaseState)
  return trailing_zeros(UInt8(x)) + 1
end

function Base.convert(::Type{DiseaseState}, x::Int64)
  return _DiseaseStateVector[x]
end

function Base.show(io::IO, x::DiseaseState)
  return print(convert(Char, x))
end

const _state_progressions = Dict{DataType, Vector{DiseaseState}}()
_state_progressions[SEIR] = [State_S; State_E; State_I; State_R]
_state_progressions[SEI] = [State_S; State_E; State_I]
_state_progressions[SIR] = [State_S; State_I; State_R]
_state_progressions[SI] = [State_S; State_I]

function advance(x::DiseaseState, ::Type{T}) where T <: EpidemicModel
  current_index = findfirst(x .== _state_progressions[T])
  return _state_progressions[T][current_index + 1]
end

function regress(x::DiseaseState, ::Type{T}) where T <: EpidemicModel
  current_index = findfirst(x .== _state_progressions[T])
  return _state_progressions[T][current_index - 1]
end

function advance!(x::DiseaseState, ::Type{T}) where T <: EpidemicModel
  x = advance(x, T)
end

function regress!(x::DiseaseState, ::Type{T}) where T <: EpidemicModel
  x = regress(x, T)
end
