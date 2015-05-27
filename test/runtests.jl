"""
runtests.jl - test Pathogen.jl functionality
Justin Angevaare
May 2015
"""
using Pathogen
using Base.Test

test_sequence=generate_sequence(700, 0.25, 0.25, 0.25, 0.25)
@test length(test_sequence) == 700

