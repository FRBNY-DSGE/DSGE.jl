using ClusterManagers, Plots, HDF5
using QuantEcon: solve_discrete_lyapunov

addprocs_sge(10,queue="background.q")
@everywhere using DSGE
srand(47)

m = SmetsWouters("ss1", testing=true)
path = dirname(@__FILE__)
filesw = "/data/dsge_data_dir/dsgejl/realtime/input_data/data"
data=readcsv("$filesw/realtime_spec=smets_wouters_hp=true_vint=110110.csv",header=true)
data = convert(Array{Float64,2}, data[1][:,2:end])
data=data'

params = h5read("$filesw/../../output_data/smets_wouters/ss0/estimate/raw/paramsmode_vint=110110.h5","params")
push!(params,m[:e_y].value,m[:e_L].value,m[:e_w].value, m[:e_π].value, m[:e_R].value, m[:e_c].value, m[:e_i].value)
update!(m,params)

system = compute_system(m)
R = system.transition.RRR
S2 = system.measurement.QQ
Φ = system.transition.TTT

N_MH = 2
n_particles = 40000
parallel = true

rand_mat = randn(size(S2,1),1)
m<=Setting(:tpf_rand_mat,rand_mat)
m<=Setting(:tpf_rstar,2.0)
m<=Setting(:tpf_N_MH, N_MH)
m<=Setting(:tpf_c,0.1)
m<=Setting(:tpf_acpt_rate,0.5)
m<=Setting(:tpf_trgt,0.25)
m<=Setting(:tpf_n_particles,n_particles)
m<=Setting(:use_parallel_workers, parallel)
m<=Setting(:x_tolerance, zero(float(0)))
m<=Setting(:tpf_deterministic, false)


s0 = zeros(size(system[:TTT])[1])
P0= nearestSPD(solve_discrete_lyapunov(Φ,R*S2*R'))
#=
Neff, lik_reg, times = tpf(m,data,system,s0,P0,0)
Neff, lik_tapered, times = tpf(m,data,system,s0,P0,1)
Neff, lik_constant, times = tpf(m,data,system,s0,P0,2)

rmprocs()

h5open("$path/../../test/reference/compare_tapering_of_error.h5","w") do file
    write(file,"lik_tapered",lik_tapered)
    write(file,"lik_reg",lik_reg)
    write(file,"lik_constant",lik_constant)
end

plotly()
plot(lik_tapered)
plot!(lik_reg)
plot!(lik_constant)
gui()
=#

Neff, lik_40000, times = tpf(m,data,system,s0,P0,0)


### BUILD NEW MODEL ###
m = SmetsWouters("ss1", testing=true)
path = dirname(@__FILE__)
filesw = "/data/dsge_data_dir/dsgejl/realtime/input_data/data"
data=readcsv("$filesw/realtime_spec=smets_wouters_hp=true_vint=110110.csv",header=true)
data = convert(Array{Float64,2}, data[1][:,2:end])
data=data'

params = h5read("$filesw/../../output_data/smets_wouters/ss0/estimate/raw/paramsmode_vint=110110.h5","params")
push!(params,m[:e_y].value,m[:e_L].value,m[:e_w].value, m[:e_π].value, m[:e_R].value, m[:e_c].value, m[:e_i].value)
update!(m,params)

system = compute_system(m)
R = system.transition.RRR
S2 = system.measurement.QQ
Φ = system.transition.TTT

N_MH = 2
n_particles = 4000
parallel = true

rand_mat = randn(size(S2,1),1)
m<=Setting(:tpf_rand_mat,rand_mat)
m<=Setting(:tpf_rstar,2.0)
m<=Setting(:tpf_N_MH, N_MH)
m<=Setting(:tpf_c,0.1)
m<=Setting(:tpf_acpt_rate,0.5)
m<=Setting(:tpf_trgt,0.25)
m<=Setting(:tpf_n_particles,n_particles)
m<=Setting(:use_parallel_workers, parallel)
m<=Setting(:x_tolerance, zero(float(0)))
m<=Setting(:tpf_deterministic, false)

s0 = zeros(size(system[:TTT])[1])
P0= nearestSPD(solve_discrete_lyapunov(Φ,R*S2*R'))
### SET UP MODEL ###

Neff, lik_4000, times = tpf(m,data,system,s0,P0,0)

rmprocs()
#=
h5open("$path/../../test/reference/compare_tapering_of_error.h5","w") do file
    write(file,"lik_tapered",lik_tapered)
    write(file,"lik_reg",lik_reg)
    write(file,"lik_constant",lik_constant)
end
=#
#=
h5open("$path/../../test/reference/compare_scaling_no_MH_change.h5","w") do file
    write(file,"lik_fifth",lik_fifth)
    write(file,"lik_tenth",lik_tenth)
end
=#

h5open("$path/../../test/reference/compare_n_particles.h5","w") do file
    write(file,"lik_4000",lik_4000)
    write(file,"lik_40000",lik_40000)
end

lik_kalman = h5read("$path/../../test/reference/varLik.h5","kalman_lik")

#=
plotly()
plot(lik_tapered)
plot!(lik_reg)
plot!(lik_constant)
gui()
=#
#=
lik_half = h5read("$path/../../test/reference/compare_scaling.h5","lik_half")
plotly()
plot(lik_fifth,label="fifth")
plot!(lik_tenth,label="tenth")
plot!(lik_half,label="half")
plot!(lik_kalman,label="kalman")
gui()
=#

plotly()
plot(lik_4000,label="4000 particles")
plot!(lik_40000,label="40000 particles")
plot!(lik_kalman,label="kalman")
gui()

nothing