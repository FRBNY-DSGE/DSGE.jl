@everywhere using DSGE
@everywhere using Distributions

path = dirname(@__FILE__())

m = Model990()
m.testing = true

ndraws = 2
sys = compute_system(m)
syses = repmat([sys],ndraws)
zend = Vector{Vector{Float64}}()
for i in 1:ndraws
    push!(zend, zeros(DSGE.n_states_augmented(m)))
end

# fcasts = DSGE.forecast(m, syses, zend)
states, observables = DSGE.forecast(m, syses, zend)

nothing
