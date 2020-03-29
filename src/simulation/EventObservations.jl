struct EventObservations{S <: DiseaseStateSequence, M <: ILM}
  infection::Vector{Float64}
  removal::Union{Nothing, Vector{Float64}}
  seq::Union{Nothing, Vector{Union{Nothing, GeneticSeq}}}
  individuals::Int64

  function EventObservations{S, M}(i::V, r::V) where {V <: Vector{Float64}, 
                                                      S <: Union{SEIR, SIR}, 
                                                      M <: TNILM}
    if length(i) != length(r)
      @error "Length of infection and removal times must be equal"
    end
    return new{S, M}(i, r, nothing, length(i))
  end

  function EventObservations{S, M}(i::V) where {
                                   V <: Vector{Float64}, 
                                   S <: Union{SEI, SI}, 
                                   M <: TNILM}
    return new{S, M}(i, nothing, nothing, length(i))
  end


  function EventObservations{S, M}(i::V, r::V, seq::VG) where {
                                   V <: Vector{Float64},
                                   VG <: Vector{Union{Nothing, GeneticSeq}},
                                   S <: Union{SEIR, SIR}, 
                                   M <: PhyloILM}
    if length(unique([length(i); length(r); length(seq)])) != 1
      @error "Length of infection times, removal times, and genetic sequence vectors must be equal"
    end
    return new{S, M}(i, r, seq, length(i))
  end

  function EventObservations{S, M}(i::V, seq::VG) where {
                                   V <: Vector{Float64},
                                   VG <: Vector{Union{Nothing, GeneticSeq}},
                                   S <: Union{SEI, SI}, 
                                   M <: PhyloILM}
    if length(i) != length(seq)
      @error "Length of infection time and genetic sequence vectors must be equal"
    end
    return new{S, M}(i, nothing, seq, length(i))
  end
end

function EventObservations{S, M}(i::Array{Float64, 2}) where {
                                 S <: Union{SEI, SI}, 
                                 M <: TNILM}
  if size(i, 2) != 1
    @error "Invalid Array dimensions for observations of a $S model"
  end
  return EventObservations{S, M}(i[:,1])
end

function EventObservations{S, M}(ir::Array{Float64, 2}) where {
                                 S <: Union{SEIR, SIR}, 
                                 M <: TNILM}
  if size(ir, 2) != 2
    @error "Invalid Array dimensions for observations of a $S model"
  end
  return EventObservations{S, M}(ir[:, 1], ir[:, 2])
end

function EventObservations{S, M}(i::Array{Float64, 2}, seq::VG) where {
                                 VG <: Vector{Union{Nothing, GeneticSeq}},
                                 S <: Union{SEI, SI}, 
                                 M <: PhyloILM}
  if size(i, 2) != 1
    @error "Invalid Array dimensions for observations of a $S model"
  end
  return EventObservations{S, M}(i[:,1], seq)
end

function EventObservations{S, M}(ir::Array{Float64, 2}, seq::VG) where {
                                 VG <: Vector{Union{Nothing, GeneticSeq}},
                                 S <: Union{SEIR, SIR}, 
                                 M <: PhyloILM}
  if size(ir, 2) != 2
    @error "Invalid Array dimensions for observations of a $S model"
  end
  return EventObservations{S, M}(ir[:, 1], ir[:, 2], seq)
end

function Base.show(io::IO, 
                   x::EventObservations{S, M}) where {
                   S <: DiseaseStateSequence, 
                   M <: ILM}
  return print(io, "$S $M observations (n=$(x.individuals))")
end