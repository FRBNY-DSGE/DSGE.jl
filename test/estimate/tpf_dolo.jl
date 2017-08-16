using AxisArrays
import Dolo

model = Dolo.yaml_import("test/estimate/sudden_stop.yaml")

variables = cat(1, model.symbols[:states], model.symbols[:controls]) # obsrved states

n_s = length(model.symbols[:states])
n_x = length(model.symbols[:controls])

sol0 = Dolo.time_iteration(model),
dr0 = sol0.dr

sim = Dolo.simulate(model, dr0, N=1,T=50)[Axis{:N}(1)]
data = sim[Axis{:V}(variables)]  # remove auxiliaries, exogenous, etc.

dr = dr0
p = model.calibration[:parameters]



function Φ(s::Vector, ϵ::Vector)
    s_ = s[1:n_s] # Dolo states
    x_ = s[(n_s+1):(n_s+n_x)]
    # first argument is actually ignored
    S_ = Dolo.transition(model, ϵ, s_, x_, ϵ, p)
    X_ = dr(ϵ, S_)
    return cat(1, S_, X_)
end

function Ψ(s::Vector, u::Vector)
    s + u
end

Φ(s_t1::Matrix, ε_t1::Matrix) = hcat([Φ(s_t1[:,i], ε_t1[:,i]) for i in 1:size(s_t1)[2]]...)
Ψ(s_t1::Matrix, u_t1::Matrix) = hcat([Ψ(s_t1[:,i], u_t1[:,i]) for i in 1:size(s_t1)[2]]...)

s0 = cat(1, model.calibration[:states, :controls]...)
e0 = model.calibration[:exogenous]


# Convert dolo distrib to Distriutions object
import Distributions
dist_ϵ = Distributions.MvNormal(model.exogenous.mu, model.exogenous.Sigma)


import DSGE: initialize_state_draws
N_part = 100
initialize_state_draws(s0, dist_ϵ, Φ, 100)
