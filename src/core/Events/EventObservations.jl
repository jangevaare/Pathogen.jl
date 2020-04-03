struct EventObservations{S <: DiseaseStateSequence, M <: ILM}
  infection::Vector{Float64}
  removal::Union{Nothing, Vector{Float64}}
  seq::Union{Nothing, Vector{Union{Nothing, GeneticSeq}}}
  start_time::Float64
  individuals::Int64

  function EventObservations{S, M}(i::V, 
                                   r::V; 
                                   start_time::Float64=0.0) where {
                                   V <: Vector{Float64}, 
                                   S <: Union{SEIR, SIR}, 
                                   M <: TNILM}
    if length(i) != length(r)
      throw(ErrorException("Length of infection and removal times must be equal"))
    end
    return new{S, M}(i, r, nothing, start_time, length(i))
  end

  function EventObservations{S, M}(i::V; 
                                   start_time::Float64=0.0) where {
                                   V <: Vector{Float64}, 
                                   S <: Union{SEI, SI}, 
                                   M <: TNILM}
    return new{S, M}(i, nothing, nothing, start_time, length(i))
  end


  function EventObservations{S, M}(i::V, 
                                   r::V, 
                                   seq::VG; 
                                   start_time::Float64=0.0) where {
                                   V <: Vector{Float64},
                                   G <: GeneticSeq,
                                   VG <: Vector{Union{Nothing, G}},
                                   S <: Union{SEIR, SIR}, 
                                   M <: PhyloILM}
    if length(unique([length(i); length(r); length(seq)])) != 1
      throw(ErrorException("Length of infection times, removal times, and genetic sequence vectors must be equal"))
    end
    return new{S, M}(i, r, seq, start_time, length(i))
  end

  function EventObservations{S, M}(i::V, 
                                   seq::VG; 
                                   start_time::Float64=0.0) where {
                                   V <: Vector{Float64},
                                   G <: GeneticSeq,
                                   VG <: Vector{Union{Nothing, G}},
                                   S <: Union{SEIR, SIR}, 
                                   M <: PhyloILM}
    if length(i) != length(seq)
      throw(ErrorException("Length of infection time and genetic sequence vectors must be equal"))
    end
    return new{S, M}(i, nothing, seq, start_time, length(i))
  end
end

function EventObservations{S, M}(i::Array{Float64, 2}; 
                                 start_time::Float64=0.0) where {
                                 S <: Union{SEI, SI}, 
                                 M <: TNILM}
  if size(i, 2) != 1
    throw(ErrorException("Invalid Array dimensions for observations of a $S model"))
  end
  return EventObservations{S, M}(i[:,1], start_time=start_time)
end

function EventObservations{S, M}(ir::Array{Float64, 2}; 
                                 start_time::Float64=0.0) where {
                                 S <: Union{SEIR, SIR}, 
                                 M <: TNILM}
  if size(ir, 2) != 2
    throw(ErrorException("Invalid Array dimensions for observations of a $S model"))
  end
  return EventObservations{S, M}(ir[:, 1], ir[:, 2], start_time=start_time)
end

function EventObservations{S, M}(i::Array{Float64, 2}, 
                                 seq::VG; 
                                 start_time::Float64=0.0) where {
                                 G <: GeneticSeq,
                                 VG <: Vector{Union{Nothing, G}},
                                 S <: Union{SEI, SI}, 
                                 M <: PhyloILM}
  if size(i, 2) != 1
    throw(ErrorException("Invalid Array dimensions for observations of a $S model"))
  end
  return EventObservations{S, M}(i[:,1], seq, start_time=start_time)
end

function EventObservations{S, M}(ir::Array{Float64, 2}, 
                                 seq::VG; 
                                 start_time::Float64=0.0) where {
                                 G <: GeneticSeq,
                                 VG <: Vector{Union{Nothing, G}},
                                 S <: Union{SEIR, SIR}, 
                                 M <: PhyloILM}
  if size(ir, 2) != 2
    throw(ErrorException("Invalid Array dimensions for observations of a $S model"))
  end
  return EventObservations{S, M}(ir[:, 1], ir[:, 2], seq, start_time=start_time)
end

function Base.show(io::IO, 
                   x::EventObservations{S, M}) where {
                   S <: DiseaseStateSequence, 
                   M <: ILM}
  return print(io, "$S $M observations (n=$(x.individuals))")
end