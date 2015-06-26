using Pathogen

phytree=Tree([generate_seq(200, 0.25, 0.25, 0.25, 0.25), generate_seq(200, 0.25, 0.25, 0.25, 0.25), generate_seq(200, 0.25, 0.25, 0.25, 0.25)],
              Vector{Bool}[[false, false], [false, true], [true]],
              Vector{Float64}[[1.2, 0.3], [1.2, 2.1], [0.8]])

treedistance(1, 2, phytree)
