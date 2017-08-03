using ClusterManagers, Plots, HDF5
using QuantEcon: solve_discrete_lyapunov
using DSGE

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
n_particles = 4000
parallel = true

rand_mat = randn(size(S2,1) ,1)
m<=Setting(:tpf_rand_mat, rand_mat)
m<=Setting(:tpf_rstar, 2.0)
m<=Setting(:tpf_N_MH, N_MH)
m<=Setting(:tpf_c, 0.1)
m<=Setting(:tpf_acpt_rate, 0.5)
m<=Setting(:tpf_target, 0.4)
m<=Setting(:tpf_n_particles, n_particles)
m<=Setting(:use_parallel_workers, parallel)
m<=Setting(:x_tolerance, zero(float(0)))
m<=Setting(:tpf_deterministic, false)

s0 = zeros(size(system[:TTT])[1])

P0 = nearestSPD(solve_discrete_lyapunov(Φ,R*S2*R'))
neff, lik_fixed = tpf(m, data, system, s0, P0)

kalman_out = DSGE.filter(m,data,system,s0,P0)
lik_kalman_no_NaN = kalman_out[:L_vec]

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

#=
Neff, lik_40000, times = tpf(m,data,system,s0,P0)

h5open("$path/../../test/reference/compare_n_particles.h5","w") do file
    write(file,"lik_40000",lik_40000)
end

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
m<=Setting(:tpf_target,0.25)
m<=Setting(:tpf_n_particles,n_particles)
m<=Setting(:use_parallel_workers, parallel)
m<=Setting(:x_tolerance, zero(float(0)))
m<=Setting(:tpf_deterministic, false)

s0 = zeros(size(system[:TTT])[1])
P0= nearestSPD(solve_discrete_lyapunov(Φ,R*S2*R'))
### SET UP MODEL ###

Neff, lik_4000, times = tpf(m,data,system,s0,P0)
=#

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
#=
h5open("$path/../../test/reference/compare_n_particles.h5","w") do file
    write(file,"lik_4000",lik_4000)
    write(file,"lik_40000",lik_40000)
end
=#
#=
h5open("$path/../../test/reference/kalman_different_errors.h5","w") do file
    write(file, "lik_kalman_tenth",lik_kalman_tenth)
    write(file, "lik_kalman_half",lik_kalman_half)
end
=#
h5open("/data/dsge_data_dir/dsgejl/interns2017/standard_vals.h5","w") do file
    write(file, "lik_kalman", lik_kalman_no_NaN)
    write(file, "lik_tpf", lik_adaptive)
end 

#lik_4000 = h5read("$path/../../test/reference/compare_n_particles.h5","lik_4000")
#lik_40000 = h5read("$path/../../test/reference/compare_n_particles.h5","lik_40000")
#lik_kalman = h5read("$path/../../test/reference/varLik.h5","kalman_lik")

#=
plotly()
plot(lik_tapered)
plot!(lik_reg)
plot!(lik_constant)
gui()
=#

#lik_half = h5read("$path/../../test/reference/compare_scaling.h5","lik_half")
#lik_fifth = h5read("$path/../../test/reference/compare_scaling_no_MH_change.h5","lik_fifth")
#lik_tenth = h5read("$path/../../test/reference/compare_scaling_no_MH_change.h5","lik_tenth")
#lik_kalman_error = h5read("$path/../../test/reference/varLik_error.h5","lik_kalman_error")
plotly()
plot(lik_adaptive,label="TPF")
plot!(lik_kalman_no_NaN,label="Kalman")
gui()

#=
plotly()
plot(lik_4000,label="4000 particles")
plot!(lik_40000,label="40000 particles")
#plot!(lik_kalman,label="kalman")
plot!(lik_kalman_error,label="kalman with error")
gui()
=#
nothing