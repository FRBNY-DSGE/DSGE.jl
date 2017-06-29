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
#=
sys=compute_system(m)
Φ=sys.transition.TTT
R=sys.transition.RRR
Q=sys.measurement.QQ
=#

###For comparing piece-by-piece with MatLab, get the matrices from the Matlab code
#(copied and pasted)

A = [0.51; 3.16; 5.45]
B = [1 0 0 -1 0 1 0 0;
     0 4 0 0 0 0 0 0;
     0 0 4 0 0 0 0 0]
 H = [0.0135 0 0;
      0 0.0865 0;
      0 0 0.2003]
R = [0.7696    1.0000   -0.6674;
    1.5482    0.0000   -1.0523;
    0.7569    0.0000    0.4677;
   -0.0000   -0.0000    0.0000;
   -0.0000    1.0000    0.0000;
    1.0000   -0.0000    0.0000;
    0.3066    0.9800   -0.2528;
    0.7947    0.0000   -0.3987;
]
S2 = [0.0576         0         0;
      0    0.4225         0;
      0         0    0.0361]
Φ = [0.0000         0   -0.5406         0    0.9800    0.7158    0.0000    0.0000;
   -0.0000         0   -0.8524         0    0.0000    1.4399    0.0000   -0.0000;
   -0.0000         0    0.3788         0    0.0000    0.7039   -0.0000    0.0000;
    1.0000         0   -0.0000         0   -0.0000   -0.0000    0.0000    0.0000;
   -0.0000         0    0.0000         0    0.9800   -0.0000   -0.0000    0.0000;
    0.0000         0    0.0000         0   -0.0000    0.9300    0.0000    0.0000;
    0.0000         0   -0.2048         0    0.9604    0.2851   -0.0000    0.0000;
    0.0000         0   -0.3229         0    0.0000    0.7390    0.0000    0.0000]

s0 = zeros(8)
P0=nearestSPD(solve_discrete_lyapunov(Φ, R*S2*R'))

#=println(size(A))
println(size(B))
println(size(H))
println(size(R))
println(size(S2))
println(size(Φ))
=#
#println(path)
df = readtable("$path/../../../../../us.txt",separator=' ')
data = convert(Matrix{Float64},df)
data=data'
#println(data)
neff, lik = tpf(m, data, s0, P0, A, B, H, R, S2, Φ)
@show neff, lik