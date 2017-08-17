using AxisArrays
using DSGE, Distributions
using Plots
import Dolo
gr()

model = Dolo.yaml_import("sudden_stop.yaml")
variables = vcat(model.symbols[:states], model.symbols[:controls]) # obsrved states

n_s = length(model.symbols[:states])
n_x = length(model.symbols[:controls])

sol0 = Dolo.time_iteration(model)
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
    return vcat(S_, X_)
end

function Ψ(s::Vector, u::Vector)
    s[2:3] + u
end
# Because we've chosen the data to only be on the middle two states
data = data[2:3,:]

# Φ(s_t1::Matrix, ε_t1::Matrix) = hcat([Φ(s_t1[:,i], ε_t1[:,i]) for i in 1:size(s_t1)[2]]...)
# Ψ(s_t1::Matrix, u_t1::Matrix) = hcat([Ψ(s_t1[:,i], u_t1[:,i]) for i in 1:size(s_t1)[2]]...)

s0 = vcat(model.calibration[:states, :controls]...)
e0 = model.calibration[:exogenous]

# Convert dolo distrib to Distriutions object
F_ϵ = Distributions.MvNormal(model.exogenous.mu, model.exogenous.Sigma)
F_u = Distributions.MvNormal(zeros(2), 1e-8*eye(2))

# Instantiating a model purely for obtaining settings for tpf algorithm
m = AnSchorfheide()

N_MH = 1
n_parts = 4000
parallel = false

m <= Setting(:tpf_n_mh_simulations, N_MH)
m <= Setting(:tpf_c_star, 0.1)
m <= Setting(:tpf_n_particles, n_parts)
m <= Setting(:use_parallel_workers, parallel)
m <= Setting(:x_tolerance, 0.)
m <= Setting(:tpf_deterministic, false)

s_init = initialize_state_draws(s0, F_ϵ, Φ, n_parts)
Neff, lik_tpf = tpf(m, data.data, Φ, Ψ, F_ϵ, F_u, s_init)
