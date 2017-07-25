using ClusterManagers, Plots, HDF5
using QuantEcon: solve_discrete_lyapunov
#=
addprocs_sge(10,queue="background.q")
@everywhere using DSGE
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

rand_mat = randn(size(S2,1),1)
m<=Setting(:tpf_rand_mat,rand_mat)

m<=Setting(:tpf_rstar,2.0)
m<=Setting(:tpf_N_MH, 10)
m<=Setting(:tpf_c,0.1)
m<=Setting(:tpf_acpt_rate,0.5)
m<=Setting(:tpf_trgt,0.25)
n_particles=4000
m<=Setting(:tpf_n_particles,n_particles)
m<=Setting(:use_parallel_workers, true)
m<=Setting(:x_tolerance, zero(float(0)))
m<=Setting(:tpf_deterministic, false)


s0 = zeros(size(system[:TTT])[1])
P0= nearestSPD(solve_discrete_lyapunov(Φ,R*S2*R'))

Neff,lik,times = tpf(m,data,system,s0,P0)
times10 = times
rmprocs()

h5open("$path/../../test/reference/timingAndWorkers10.h5","w") do file
    write(file,"times10",times10)
end

@show procs()
=#
addprocs_sge(30,queue="background.q")
@show procs()

@everywhere using DSGE
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

rand_mat = randn(size(S2,1),1)
m<=Setting(:tpf_rand_mat,rand_mat)

m<=Setting(:tpf_rstar,2.0)
m<=Setting(:tpf_N_MH, 10)
m<=Setting(:tpf_c,0.1)
m<=Setting(:tpf_acpt_rate,0.5)
m<=Setting(:tpf_trgt,0.25)
n_particles=4000
m<=Setting(:tpf_n_particles,n_particles)
m<=Setting(:use_parallel_workers, true)
m<=Setting(:x_tolerance, zero(float(0)))
m<=Setting(:tpf_deterministic, false)


s0 = zeros(size(system[:TTT])[1])
P0= nearestSPD(solve_discrete_lyapunov(Φ,R*S2*R'))

Neff, lik, times = tpf(m,data,system,s0,P0)
times30 = times
rmprocs()
h5open("$path/../../test/reference/timingAndWorkers30.h5","w") do file
    write(file,"times30",times30)
end
plotly()
plot(times10)
plot!(times30)
gui()