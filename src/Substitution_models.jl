function JC69(p::Tuple)
  """
  Returns the JC69 (Jukes and Cantor, 1969) substitution rate matrix with parameter μ
  """
  @assert 0 < p[1] ["Invalid parameter value"]
  [[0      p[1]/4  p[1]/4  p[1]/4]
   [p[1]/4 0       p[1]/4  p[1]/4]
   [p[1]/4 p[1]/4  0       p[1]/4]
   [p[1]/4 p[1]/4  p[1]/4  0]]
end

