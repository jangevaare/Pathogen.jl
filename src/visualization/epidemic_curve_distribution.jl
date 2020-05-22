function _epidemic_curve_distribution(x::Vector{Events{S}},
                                      state::DiseaseState,
                                      times;
                                      credibleinterval::Float64=0.95) where {
                                      S <: DiseaseStateSequence}
  !(0.0 <= credibleinterval <= 1.0) && error("`credibleinterval` must be between 0 and 1")
  ci = (1-credibleinterval)/2
  lwr  = Array{Tuple{Float64,Float64}, 1}(undef, length(times))
  upr  = Array{Tuple{Float64,Float64}, 1}(undef, length(times))
  for t in eachindex(times)
    y = _count_by_state.(x, Ref(state), times[t])
    q = quantile(y, [ci, 1-ci])
    lwr[t] = (times[t], q[1])
    upr[t] = (times[t], q[2])
  end
  return [lwr; reverse(upr)]
end



@recipe function f(x::Vector{Events{S}},
                   state::DiseaseState,
                   times;
                   credibleinterval=0.95)  where {
                   S <: DiseaseStateSequence}
  xguide --> "Time"
  yguide --> "N"
  label --> convert(Char, state)
  legend --> :topright
  fillcolor --> _state_color(state)
  linecolor --> _state_color(state)
  fillalpha --> 0.5
  linealpha --> 0.7
  seriestype := :shape
  _epidemic_curve_distribution(x, state, times, credibleinterval=credibleinterval)
end

@recipe function f(x::Vector{Events{S}},
                   times;
                   credibleinterval=0.95) where {
                   S<: DiseaseStateSequence}
  for s in convert(Tuple, S)
    @series begin
      credibleinterval := credibleinterval
      x, s, times
    end
  end
end