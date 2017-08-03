using DSGE, HDF5, DataFrames
using QuantEcon: solve_discrete_lyapunov
using Plots

m = AnSchorfheide(testing=true)
path=dirname(@__FILE__)
df = readtable("$path/../../test/reference/us.txt", header=false, separator=' ')
data = convert(Matrix{Float64}, df)
data = data'

params = [2.09,0.98,2.25,0.65,0.34,3.16,0.51,0.81,0.98,0.93,0.19,0.65,0.24,0.12,0.29,0.45]
update!(m,params)

system = compute_system(m)
R = system.transition.RRR
S2 = system.measurement.QQ
phi = system.transition.TTT

m<=Setting(:tpf_rstar, 2.0)
m<=Setting(:tpf_N_MH, 2)
m<=Setting(:tpf_c, 0.1)
m<=Setting(:tpf_acpt_rate,0.5)
m<=Setting(:tpf_trgt,0.25)
n_particles=4000
m<=Setting(:tpf_n_particles,n_particles)
m<=Setting(:use_parallel_workers, false)
m<=Setting(:x_tolerance, zero(float(0)))
m<=Setting(:tpf_deterministic,false)

s0=zeros(size(system[:TTT])[1])
P0=nearestSPD(solve_discrete_lyapunov(phi,R*S2*R'))

kalman_out = DSGE.filter(m,data,s0,P0)
kalman = kalman_out[:L_vec]

tpf_liks = zeros(80,15)
for i=1:15
    Neff, lik, times = tpf(m,data,system,s0,P0,0)
    tpf_liks[:,i] = lik
end

plotly()
plot(kalman, label="kalman")
for i=1:15
    plot!(tpf_liks[:,i], label="tpf $(i)")
end
gui()


