"""
Return the transmission pathway leading to an individual
"""
function pathwayto(individual::Int64,
                   events::Events,
                   depth=Inf::Real)
  if isnan(events.exposed[individual])
    path = Int64[]
  else
    path = [individual]
    while length(path) <= depth && sum(events.network[2][:, path[end]]) == 1
      push!(path, findfirst(events.network[2][:, path[end]]))
    end
  end
  return path
end


"""
Return all transmission pathways leading to specified individuals
"""
function pathwayto(individuals::Vector{Int64},
                   events::Events,
                   depth=Inf::Real)
  paths = Array[Int64[]]
  for i in individuals
    push!(paths, pathwayto(individual, events, depth))
  end
  return paths
end


"""
Return all transmission pathways leading to individuals
"""
function pathwayto(events::Events,
                   depth=Inf::Real)
  return pathwayto(1:length(events.exposed), events, depth)
end


"""
Return the transmission pathway leading from an individual
"""
function pathwayfrom(individual::Int64,
                     events::Events,
                     depth=Inf::Int64)
  if isnan(events.exposed[individual])
   path = Int64[]
  else
   path = [individual]
   pathlengths = [0]
   while depth >= length(path) > pathlengths[end]
     push!(pathlengths, length(path))
     for j in path[(pathlengths[end-1]+1):pathlengths[end]]
       append!(path, find(events.network[2][j, :]))
     end
   end
  end
  return path
end


"""
Return all transmission pathways leading from specified individuals
"""
function pathwayfrom(individuals::Vector{Int64},
                     events::Events,
                     depth=Inf::Int64)
  paths = Array[Int64[]]
  for i in individuals
    push!(paths, pathwayfrom(i))
  end
  return paths
end


"""
Return all transmission pathways leading from individuals
"""
function pathwayfrom(events::Events,
                      depth=Inf::Int64)
  return pathwayfrom(1:length(events.exposed), events, depth)
end
