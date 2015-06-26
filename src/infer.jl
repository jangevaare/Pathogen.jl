"""
infer.jl - pathogen evolution and transmission dynamic inference tools
Justin Angevaare
June 2015
"""

function branchloglikelihood(seq1::Nucleotide2bitSeq, seq2::Nucleotide2bitSeq, branchdistance::Float64, substitution_matrix::Array)
  """
  Log likelihood for any two aligned sequences, a specified distance apart on a phylogenetic tree
  """
  @assert(length(seq1) == length(seq2), "Sequences not aligned")
  ll = 0
  for i = 1:length(seq1)
    base1 = convert(Int64, seq1[i])
    for base2 = 1:4
      if base2 == convert(Int64, seq2[i])
        if base1 != base2
          ll += log(1 - exp(substitution_matrix[base1, base2] .* branchdistance))
        end
      else
        ll += substitution_matrix[base1, base2] .* branchdistance
      end
    end
  end
  return ll
end


