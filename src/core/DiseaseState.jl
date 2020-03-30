primitive type DiseaseState 8 end

DiseaseState(x::UInt8) = reinterpret(DiseaseState, x)
UInt8(x::DiseaseState) = reinterpret(UInt8, x)

# `DieaseState` const representations
const State_S = DiseaseState(0b0001)
const State_E = DiseaseState(0b0010)
const State_I = DiseaseState(0b0100)
const State_R = DiseaseState(0b1000)

const DiseaseStates = Vector{DiseaseState}

const _DiseaseStateVector = [State_S; State_E; State_I; State_R]
const _DiseaseStateCharVector = ['S'; 'E'; 'I'; 'R']

# `DiseaseState` conversion to/from `Char`
function Base.convert(::Type{Char}, x::DiseaseState)
  return _DiseaseStateCharVector[_DiseaseStateVector .== Ref(x)][1]
end

function Base.show(io::IO, x::DiseaseState)
  return print(io, convert(Char, x))
end

function Base.convert(::Type{DiseaseStates}, x::Type{SEIR})
  return [State_S, State_E, State_I, State_R]
end

function Base.convert(::Type{DiseaseStates}, x::Type{SEI})
  return [State_S, State_E, State_I]
end

function Base.convert(::Type{DiseaseStates}, x::Type{SIR})
  return [State_S, State_I, State_R]
end

function Base.convert(::Type{DiseaseStates}, x::Type{SI})
  return [State_S, State_I]
end



# function regress(x::DiseaseState, ::Type{S}) where {S <: DiseaseStateSequence}
#   states = convert(DiseaseStates, S)
#   current_index = findfirst(Ref(x) .== states)
#   return states[current_index - 1]
# end

# function Base.length(x::DiseaseState)
#   return 1
# end

# function Base.convert(::Type{DiseaseState}, x::Char)
#   return _DiseaseStateVector[_DiseaseStateCharVector .== Ref(x)][1]
# end

# # `DiseaseState` conversion to/from `BitArray`
# function Base.convert(::Type{BitArray{1}}, x::DiseaseState)
#   return BitArray(_DiseaseStateVector .== Ref(x))
# end

# function Base.convert(::Type{DiseaseState}, x::BitArray{1})
#   return _DiseaseStateVector[x]
# end

# # `DiseaseState` conversion to/from `Int64`
# function Base.convert(::Type{Int64}, x::DiseaseState)
#   return trailing_zeros(UInt8(x)) + 1
# end

# function Base.convert(::Type{DiseaseState}, x::Int64)
#   return _DiseaseStateVector[x]
# end


# function advance!(x::DiseaseState, ::Type{S}) where {S <: DiseaseStateSequence}
#   x = advance(x, M)
# end

# function regress!(x::DiseaseState, ::Type{S}) where {S <: DiseaseStateSequence}
#   x = regress(x, M)
# end