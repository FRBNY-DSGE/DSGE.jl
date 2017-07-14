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
	rand_mat = randn(3,1)
    	
	h5open("$path/../reference/mutationRandomMatrix.h5","w") do file
            write(file,"rand_mat",rand_mat)
    	end
        
        # Set variables within system
    	transition_equation = Transition(Φ, R)
    	measurement_equation = Measurement(B,squeeze(A,2),S2,H,rand_mat,R)
    	system = System(transition_equation, measurement_equation)
        
        # params = [2.09, 0.98, 2.25, 0.65, 0.34, 3.16, 0.51, 0.81, 0.98, 0.93, 0.19, 0.65, 0.24,0.12,0.29,0.45]
        # update!(m,params)
    else
        # If not testing, compute system in Julia, get better starting parameters s.t. code runs
        file = "$path/../reference/optimize.h5"

        x0 = h5read(file,"params")
        data = h5read(file, "data")'
        
        #filesw = "/data/dsge_data_dir/dsgejl/realtime/input_data/data"
        #data = readcsv("$filesw/realtime_spec=smets_wouters_hp=true_vint=110110.csv",header=true)
        #data = convert(Array{Float64,2}, data[1][:,2:end])
        #data=data'
       
        minimizer = h5read(file,"minimizer")
        update!(m,x0)
        x0=Float64[p.value for p in m.parameters]
#        params = h5read("$filesw/../../output_data/smets_wouters/ss0/estimate/raw/paramsmode_vint=110110.h5","params")
        out, H = optimize!(m, data; iterations=200)
        params = out.minimizer
        update!(m,params)

        system = compute_system(m)
        R = system.transition.RRR
        S2 = system.measurement.QQ
        Φ = system.transition.TTT

    end
    return system, data, Φ, R, S2
end


# Set up model

# An Schorfheide model
custom_settings = Dict{Symbol, Setting}(
    :date_forecast_start => Setting(:date_forecast_start, quartertodate("2015-Q4")))
m = AnSchorfheide(custom_settings = custom_settings, testing = true)

# Smets Wouters model
#custom_settings = Dict{Symbol, Setting}(
#    :date_forecast_start => Setting(:date_forecast_start, quartertodate("2011-Q1")))
#m = SmetsWouters(custom_settings = custom_settings, testing = true)


path=dirname(@__FILE__)

# For comparison test
good_likelihoods = h5open("$path/../reference/tpf_test_likelihoods.h5","r") do file
    read(file, "test_likelihoods")
end

# Tuning Parameters
m<=Setting(:tpf_rstar,2.0)
m<=Setting(:tpf_c,0.1)
m<=Setting(:tpf_acpt_rate,0.5)
m<=Setting(:tpf_trgt,0.25)
m<=Setting(:tpf_n_mh_simulations,2)
deterministic = false
m<=Setting(:tpf_deterministic,deterministic)
m<=Setting(:n_presample_periods,2)

# Parallelize
m<=Setting(:use_parallel_workers,true)

# Set tolerance in fzero
m<=Setting(:tpf_x_tolerance,1e-3)
#m<=Setting(:tpf_x_tolerance, zero(float(0)))

# Set number of particles
#n_particles = 500
#m<=Setting(:tpf_n_particles, n_particles)

#sys, data, Φ, R, S2  = setup(false)

#s0 = zeros(8)
#P0 = nearestSPD(solve_discrete_lyapunov(Φ, R*S2*R'))
#m<=Setting(:tpf_deterministic, true)
#tic()
#neff, lik = tpf(m, data,sys, s0, P0)
#toc()

# Test 4000 particles, testing = true
n_particles = 500
deterministic = false
m<=Setting(:tpf_n_particles, n_particles)

sys, data, Φ, R, S2 = setup(deterministic)
s0 = zeros(size(sys[:TTT])[1])
P0 = nearestSPD(solve_discrete_lyapunov(Φ, R*S2*R'))
m<=Setting(:tpf_deterministic, deterministic)

tic()
neff, lik = tpf(m, data, sys, s0, P0)
toc()

if (n_particles == 4000) & deterministic
    @test good_likelihoods == lik
    @show good_likelihoods
    @show lik
    println("Test passed for 4000 particles in testing mode.")
end



#####The following code regenerates the test comparison that we use to compare. DO NOT RUN (unless you are sure that the new tpf.jl is correct).
# Seeded, deterministic resampling; fixed tempering schedule of 0.25->0.5->1
# h5open("$path/../reference/tpf_test_likelihoods.h5","w") do file
#     write(file,"test_likelihoods",lik)
# end