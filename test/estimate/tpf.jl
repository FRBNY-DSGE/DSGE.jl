using DSGE, HDF5, DataFrames, ClusterManagers, Plots
using QuantEcon: solve_discrete_lyapunov

# Establish path
path = dirname(@__FILE__)

# Set up model
function setup_model(model_type::String, n_particles::Int64, adaptive::Bool, parallel::Bool)

    if (model_type=="AnSchorfheide")
        # An Schorfheide model
        custom_settings = Dict{Symbol, Setting}(:date_forecast_start => 
                          Setting(:date_forecast_start, quartertodate("2015-Q4")))
        m = AnSchorfheide(custom_settings = custom_settings, testing = true)
    else
        # Smets Wouters model
        custom_settings = Dict{Symbol, Setting}(:date_forecast_start => 
                          Setting(:date_forecast_start, quartertodate("2011-Q2")))
        m = SmetsWouters("ss1",custom_settings = custom_settings, testing = true)
        m<= Setting(:date_conditional_end, quartertodate("2011-Q1"))
        m<= Setting(:date_forecast_start, quartertodate("2011-Q2"))
    end

    # Tuning Parameters
    m<=Setting(:tpf_r_star, 2.0)
    m<=Setting(:tpf_c_star, 0.1)
    m<=Setting(:tpf_accept_rate,0.5)
    m<=Setting(:tpf_target, 0.25)
    m<=Setting(:tpf_n_mh_simulation, 2)
    m<=Setting(:n_presample_periods, 2)
    m<=Setting(:tpf_adaptive, adaptive)
    m<=Setting(:use_parallel_workers, parallel)
    m<=Setting(:tpf_x_tolerance, zero(float(0)))
    m<=Setting(:tpf_n_particles, n_particles)
    return m
end

# Function sets up model based on type and φ schedule
function optimize_setup(m::AbstractModel, adaptive::Bool)
    
    if typeof(m) == AnSchorfheide{Float64} 
        if !adaptive
            println("AnSchorfheide Model, fixed φ schedule.")
            # Seed random number generator
            srand(1234)
            # Load in us.txt data from schorfheide
            df = readtable("$path/../reference/us.txt", header=false, separator=' ')
            data = convert(Matrix{Float64}, df)
            data = data'

            # Read in matrices from Schorfheide Matlab code
            file = "$path/../reference/matlab_variable_for_testing.h5"
            A  = h5read(file, "A")
            B  = h5read(file, "B")
            H  = h5read(file, "H")
            R  = h5read(file, "R")
            S2 = h5read(file, "S2")
            Φ  = h5read(file, "Phi")
            
            # Set variables within system
    	    transition_equation = Transition(Φ, R)
    	    measurement_equation = Measurement(B,squeeze(A,2),S2,H,R,R)
    	    system = System(transition_equation, measurement_equation)
            
            # Parameters given in Schorfheide's MATLAB code
            params = [2.09, 0.98, 2.25, 0.65, 0.34, 3.16, 0.51, 0.81, 0.98, 
                       0.93, 0.19, 0.65, 0.24,0.12,0.29,0.45]
            update!(m,params)
        else
            println("AnSchorfheide model, adaptive φ schedule.")
            # Load in us.txt data from Schorfheide
            df = readtable("$path/../reference/us.txt", header=false, separator=' ')
            data = convert(Matrix{Float64},df)
            data = data'
 
            # Parameters given in Schorfheide's MATLAB code
            params = [2.09, 0.98, 2.25, 0.65, 0.34, 3.16, 0.51, 0.81, 0.98, 
                       0.93, 0.19, 0.65, 0.24,0.12,0.29,0.45]
            update!(m,params)
        end
    else
        # Smets Wouters
        println("SmetsWouters, adaptive schedule is $adaptive.")
        filesw = "/data/dsge_data_dir/dsgejl/realtime/input_data/data"
        data = readcsv("$filesw/realtime_spec=smets_wouters_hp=true_vint=110110.csv",header=true)
        data = convert(Array{Float64,2}, data[1][:,2:end])
        data = data'
    
        params = h5read("$filesw/../../output_data/smets_wouters/ss0/estimate/raw/paramsmode_vint=110110.h5","params")
        push!(params, m[:e_y].value, m[:e_L].value, m[:e_w].value, m[:e_π].value, 
                          m[:e_R].value, m[:e_c].value, m[:e_i].value)     
        update!(m,params)
    end 
        
    system = compute_system(m)
    R = system[:RRR]
    S2 = system[:QQ]
    Φ = system[:TTT]

    return m, system, data, Φ, R, S2
end

# Set parameters for testing
adaptive = true
parallel = true
n_particles = 4000

if parallel
    my_procs = addprocs_sge(10,queue="background.q")
    @everywhere using DSGE
end

### Tests:
m = setup_model("AnSchorfheide", n_particles, adaptive, parallel)
m, sys, data, Φ, R, S2 = optimize_setup(m, adaptive)
s0 = zeros(size(sys[:TTT])[1])
P0 = nearestSPD(solve_discrete_lyapunov(Φ, R*S2*R'))

tic()
neff, lik = tpf(m, data, sys, s0, P0)
total_time = toc()

# h5open("$path/../reference/output_likelihoods_ansch.h5","w") do f
#     write(f,"julia_likelihoods",lik)
# end

type_m = typeof(m)
N_MH = get_setting(m, :tpf_n_mh_simulations)
println("$n_particles Particles, $type_m, N_MH = $N_MH, elapsed time: $total_time seconds\n")

# For comparison test
good_likelihoods_adaptive = h5read("$path/../reference/tpf_test_likelihoods.h5","test_likelihoods")
good_likelihoods_random = h5read("$path/../reference/tpf_test_likelihoods_random.h5", "test_likelihoods")

if (n_particles==4000) & (adaptive) & (typeof(m) == AnSchorfheide{Float64})
    @show lik 
    @show good_likelihoods_random
    @test_matrix_approx_eq lik good_likelihoods_random
    println("Test passed for AnSchorfheide with 4000 particles and adaptive φ schedule.")
end

if ((n_particles == 4000) & !adaptive) & (typeof(m) == AnSchorfheide{Float64})
    @test lik == good_likelihoods
    println("Test passed for AnSchorfheide with 4000 particles and fixed φ schedule.")
end

#####The following code regenerates the test comparison that we use to compare. DO NOT RUN (unless you are sure that the new tpf.jl is correct).
#Seeded, adaptive resampling; fixed tempering schedule of 0.25->0.5->1
# h5open("$path/../reference/tpf_test_likelihoods_random.h5","w") do file
#     write(file,"test_likelihoods",lik)
# end

rmprocs(my_procs)

