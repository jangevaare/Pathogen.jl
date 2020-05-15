function _epidemic_curve_distribution(x::Vector{Events{S}},
                                      state::DiseaseState,
                                      times;
                                      credibleinterval::Float64=0.95) where {
                                      S <: DiseaseStateSequence}
  !(0.0 <= credibleinterval <= 1.0) && error("`credibleinterval` must be between 0 and 1")
  ci = (1-credibleinterval)/2
  lwr  = Array{Float64}(undef, length(times))
  upr  = Array{Float64}(undef, length(times))
  mn = Array{Float64}(undef, length(times))
  for t in eachindex(times)
    y = _count_by_state.(x, Ref(state), times[t])
    lwr[t], upr[t] = Tuple(quantile(y, [ci, 1-ci]))
    mn[t] = mean(y)
  end
  return mn, lwr, upr
end

@recipe function f(x::Vector{Events{S}},
                   state::DiseaseState,
                   times;
                   credibleinterval=0.95)  where {
                   S<: DiseaseStateSequence}
  mn, lwr, upr = _epidemic_curve_distribution(x, state, times, credibleinterval=credibleinterval)
  xguide --> "Time"
  yguide --> "N"
  label --> convert(Char, state)
  legend --> :none
  fillalpha --> 0.5
  ribbon := (mn - lwr, upr - mn)
  times, mn
end

@recipe function f(x::Vector{Events{S}},
                   times;
                   credibleinterval=0.95) where {
                   S<: DiseaseStateSequence}
  color --> _state_colors(S)
  for s in convert(Tuple, S)
    @series begin
      credibleinterval := credibleinterval
      x, s, times
    end
  end
end