using Pathogen
using Base.Test

nucleotide_convert("AGCU") == [1, 2, 3, 4]
nucleotide_revert([1,2,3,4]) == ["A", "G", "C", "U"]

sum(JC69((3,))) == 9
