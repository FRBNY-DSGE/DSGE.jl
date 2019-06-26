
using BenchmarkTools
using ClusterManagers
using HDF5, JLD, Base.Test
using DataFrames
using DSGE
import DSGE.update!
using QuantEcon: solve_discrete_lyapunov, solve_discrete_riccati

include("test.jl")
include("kalman_filter.jl")
include("ar1_ct_kalmanfunction.jl")

data, T, R, C, Q, Z, D, E, n_states = get_oneassethank_statespace()

# Generation of the initial state draws
s_0 = zeros(n_states)
P_0 = solve_discrete_lyapunov(T, R*Q*R')

##

Φ(s_t::Vector{Float64}, ϵ_t::Vector{Float64}) = TTT*s_t + RRR*ϵ_t + CCC
Ψ(s_t::Vector{Float64}, u_t::Vector{Float64}) = ZZ*s_t + DD + u_t

F_ϵ = Distributions.MvNormal(zeros(size(QQ, 1)), QQ)
F_u = Distributions.MvNormal(zeros(size(HH, 1)), HH)

# Tuning of the tempered particle filter algorithm

tuning = Dict(:r_star => 2., :c => 0.3, :accept_rate => 0.4, :target => 0.4,
              :xtol => 0., :resampling_method => :systematic, :N_MH => 1,
              :n_particles => 1000, :n_presample_periods => 0,
              :adaptive => true, :allout => true, :parallel => false)

# Generation of the initial state draws

n_states = n_states_augmented(m)
s0 = zeros(n_states)
P0 = solve_discrete_lyapunov(TTT, RRR*QQ*RRR')
U, E, V = svd(P0)
s_init = s0 .+ U*diagm(sqrt.(E))*randn(n_states, tuning[:n_particles])
tempered_particle_filter(data, Φ, Ψ, F_ϵ, F_u, s_init; tuning...)

##

#Kalman
norm_P_T, ch_ll, truelik = compute_values(data, T, R, C, Q, Z, D, E, s_0, P_0)

kaloutput = Dict{Symbol, Any}()

kaloutput[:norm_P_T] = norm_P_T
kaloutput[:ch_ll] = ch_ll
kaloutput[:truelik] = truelik