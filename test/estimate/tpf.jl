using DSGE, HDF5, DataFrames, ClusterManagers, Plots
using QuantEcon: solve_discrete_lyapunov

function setup(testing::Bool)
    if testing
        srand(1234)
        df = readtable("$path/../reference/us.txt",header=false, separator=' ')
        data = convert(Matrix{Float64},df)
        data=data'
        file = "$path/../reference/matlab_variable_for_testing.h5"
        A = h5read(file, "A")
        B = h5read(file, "B")
        H = h5read(file, "H")
        R = h5read(file, "R")
        S2 = h5read(file,"S2")
        Φ = h5read(file, "Phi")
        
        m<=Setting(:DD,A[:,1])
        m<=Setting(:ZZ,B)
        m<=Setting(:RRR,R)
        m<=Setting(:TTT,Φ)
        m<=Setting(:EE,H)
        m<=Setting(:tpf_S2,S2)
        C = zeros(8,1)
        m<=Setting(:CCC, C)
 #       params = [2.09, 0.98, 2.25, 0.65, 0.34, 3.16, 0.51, 0.81, 0.98, 0.93, 0.19, 0.65, 0.24,0.12,0.29,0.45]
      #  update!(m,params)
        return data, Φ, R, S2
    else
        #If not testing, compute system in Julia, get better starting parameters s.t. code runs
        sys=compute_system(m)
        R = sys.transition.RRR
        S2 = sys.measurement.QQ
        Φ = sys.transition.TTT

        file = "$path/../reference/optimize.h5"

        x0 = h5read(file,"params")
        data = h5read(file, "data")'

        minimizer = h5read(file,"minimizer")
        update!(m,x0)

        x0=Float64[p.value for p in m.parameters]
        out, H = optimize!(m, data; iterations=500)
        params = out.minimizer
        @show params
        update!(m,params)
        return data, Φ, R, S2
    end
end


# Set up model
custom_settings = Dict{Symbol, Setting}(
    :date_forecast_start => Setting(:date_forecast_start, quartertodate("2015-Q4")))
m = AnSchorfheide(custom_settings = custom_settings, testing = true)

path=dirname(@__FILE__)

# For comparison test
good_likelihoods = h5open("$path/../reference/tpf_test_likelihoods.h5","r") do file
    read(file, "test_likelihoods")
end

m<=Setting(:tpf_rstar,2.0)
m<=Setting(:tpf_c,0.1)
m<=Setting(:tpf_acpt_rate,0.5)
m<=Setting(:tpf_trgt,0.25)
m<=Setting(:tpf_N_MH,2)
n_particles=500
m<=Setting(:tpf_n_particles, n_particles)

# Parallelize
m<=Setting(:use_parallel_workers,true)

# Set tolerance in fzero
#m<=Setting(:x_tolerance,1e-3)
m<=Setting(:x_tolerance, zero(float(0)))
# Input parameters

# data,Φ, R, S2  = setup(false)
# s0 = zeros(8)
# P0 = nearestSPD(solve_discrete_lyapunov(Φ, R*S2*R'))
# neff, lik = tpf(m, data, s0, P0, testing=false)

data, Φ, R, S2 = setup(true)
s0 = zeros(8)
P0 = nearestSPD(solve_discrete_lyapunov(Φ, R*S2*R'))
neff, lik = tpf(m, data, s0, P0, testing=true)
n_particles=4000
m<=Setting(:tpf_n_particles, n_particles)

@test good_likelihoods == lik



#####The following code regenerates the test comparison that we use to compare. DO NOT RUN (unless you are sure that the new tpf.jl is correct).
# Seeded, deterministic resampling; fixed tempering schedule of 0.25->0.5->1
# h5open("$path/../reference/tpf_test_likelihoods.h5","w") do file
#     write(file,"test_likelihoods",lik)
# end



#I believe the following code is extraneous
#Test if the code runs without error
#=s0 = zeros(8)
ε0 = zeros(3)=#
#data = h5open("$path/../reference/mutation_RWMH.h5","r") do file
#    read(file, "data")
#end

#=
sys=compute_system(m)
Φ=sys.transition.TTT
R=sys.transition.RRR
Q=sys.measurement.QQ
=#


