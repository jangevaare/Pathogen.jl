"""
runtests.jl - test Pathogen.jl functionality
Justin Angevaare
May 2015
"""
using Pathogen
using Base.Test

@test nucleotide_convert("AGCU") == [1, 2, 3, 4]
@test nucleotide_revert([1,2,3,4]) == ["A", "G", "C", "U"]

test_sequence=generate_sequence(700, 0.25, 0.25, 0.25, 0.25)
@test length(test_sequence) == 700
@test sum(test_sequence .!= mutate_sequence(test_sequence, JC69((1,),100))) > 0
