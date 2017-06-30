using DSGE, HDF5, DataFrames
using QuantEcon: solve_discrete_lyapunov

m = AnSchorfheide(testing=true)
m<=Setting(:data_forecast_start,quartertodate("2015-Q4"))
m<=Setting(:tpf_rstar,2.0)
m<=Setting(:tpf_c,0.1)
m<=Setting(:tpf_acpt_rate,0.5)
m<=Setting(:tpf_trgt,0.25)
m<=Setting(:tpf_N_MH,2)
m<=Setting(:tpf_num_particles,100)

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

A = h5open("$path/../reference/matlab_variable_for_testing.h5","r") do file
    read(file,"A")
end

B = h5open("$path/../reference/matlab_variable_for_testing.h5","r") do file
    read(file,"B")
end

H = h5open("$path/../reference/matlab_variable_for_testing.h5","r") do file
    read(file,"H")
end

R = h5open("$path/../reference/matlab_variable_for_testing.h5","r") do file
    read(file,"R")
end

S2 = h5open("$path/../reference/matlab_variable_for_testing.h5","r") do file
    read(file,"S2")
end

Φ = h5open("$path/../reference/matlab_variable_for_testing.h5","r") do file
    read(file,"Phi")
end


###For comparing piece-by-piece with MatLab, get the matrices from the Matlab code
#(copied and pasted)


s0 = zeros(8)
P0=nearestSPD(solve_discrete_lyapunov(Φ, R*S2*R'))
df = readtable("$path/../../../../../us.txt",header=false, separator=' ')
data = convert(Matrix{Float64},df)
data=data'
neff, lik, sVec= tpf(m, data, s0, P0, A, B, H, R, S2, Φ)

h5open("$path/../reference/output_likelihoods.h5","w") do file
    write(file,"julia_likelihoods",lik)
end
@show sVec