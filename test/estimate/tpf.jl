using DSGE, HDF5, DataFrames
using QuantEcon: solve_discrete_lyapunov

m = AnSchorfheide(testing=true)
m<=Setting(:data_forecast_start,quartertodate("2015-Q4"))
m<=Setting(:tpf_rstar,3.0)
m<=Setting(:tpf_c,0.1)
m<=Setting(:tpf_acpt_rate,0.5)
m<=Setting(:tpf_trgt,0.25)
m<=Setting(:tpf_N_MH,3)
m<=Setting(:tpf_numParticles,100)

srand(1234)

path=dirname(@__FILE__)

#Test if the code runs without error
#=s0 = zeros(8)
ε0 = zeros(3)
data = h5open("$path/../reference/mutation_RWMH.h5","r") do file
    read(file, "data")
end
=#

sys=compute_system(m)
Φ=sys.transition.TTT
R=sys.transition.RRR
Q=sys.measurement.QQ

s0 = zeros(8)
P0=nearestSPD(solve_discrete_lyapunov(Φ, R*Q*R'))

println(path)
df = readtable("$path/../../../../../us.txt",separator=' ')
data = convert(Matrix{Float64},df)
data=data'
println(data)
neff, lik = tpf(m, data, s0, P0)