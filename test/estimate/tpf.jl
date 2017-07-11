using DSGE, HDF5, DataFrames, ClusterManagers, Plots
using QuantEcon: solve_discrete_lyapunov

#Set up model
custom_settings = Dict{Symbol, Setting}(
    :date_forecast_start => Setting(:date_forecast_start, quartertodate("2015-Q4")))
m = AnSchorfheide(custom_settings = custom_settings, testing = true)

path=dirname(@__FILE__)

#For comparison test
good_likelihoods = h5open("$path/../reference/tpf_test_likelihoods.h5","r") do file
    read(file, "test_likelihoods")
end

#Set settings
m<=Setting(:tpf_rstar,2.0)
m<=Setting(:tpf_c,0.1)
m<=Setting(:tpf_acpt_rate,0.5)
m<=Setting(:tpf_trgt,0.25)
m<=Setting(:tpf_N_MH,2)
numParticles=500
m<=Setting(:tpf_n_particles, numParticles)

# Parallelize
m<=Setting(:use_parallel_workers,false)

# Set testing condition
testing = false

#If testing, get us.txt data from Shorfheide package, read in Matlab matrices, seed random number generator.
if testing
    srand(1234)
    df = readtable("$path/../../../../../us.txt",header=false, separator=' ')
    data = convert(Matrix{Float64},df)
    data=data'
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
    m<=Setting(:DD,A[:,1])
    m<=Setting(:ZZ,B)
    m<=Setting(:RRR,R)
    m<=Setting(:TTT,Φ)
    m<=Setting(:EE,H)
    m<=Setting(:tpf_S2,S2)
else
    #If not testing, compute system in Julia, get better starting parameters to that the whole code runs
    sys=compute_system(m)
    R = sys.transition.RRR
    S2=sys.measurement.QQ
    Φ = sys.transition.TTT
    file = "$path/../reference/optimize.h5"
    x0 = h5read(file,"params")
    data = h5read(file, "data")'
    minimizer = h5read(file,"minimizer")
    minimum = h5read(file,"minimum")
    H_expected = h5read(file,"H")
    update!(m,x0)
    x0=Float64[p.value for p in m.parameters]
    out, H = optimize!(m,data; iterations=500)
    params = out.minimizer
    update!(m,params)
end

# Input parameters
s0 = zeros(8)
P0 = nearestSPD(solve_discrete_lyapunov(Φ, R*S2*R'))

tic()
neff, lik = tpf(m, data, s0, P0, testing=testing)
toc()

##Plotting: this will need to move, but in test file for now
plotly()
plot(lik)
gui()

#neff, lik = tpf(m, data, s0, P0, testing=0, parallel=0)
#neff, lik = tpf(m, data, s0, P0, testing=1, parallel=0)

if n_particles == 4000 & testing
    @test good_likelihoods == lik
end

h5open("$path/../reference/output_likelihoods.h5","w") do file
    write(file,"julia_likelihoods",lik)
end


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
