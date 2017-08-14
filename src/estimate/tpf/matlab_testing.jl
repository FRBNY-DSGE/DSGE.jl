using ClusterManagers, DSGE, DataFrames, Plots, HDF5
using QuantEcon: solve_discrete_lyapunov

#srand(47)
custom_settings = Dict{Symbol, Setting}(:date_forecast_start =>
                                        Setting(:date_forecast_start, quartertodate("2015-Q4")))
m = AnSchorfheide(custom_settings = custom_settings, testing=true)
path = dirname(@__FILE__)


N_MH = 2
n_particles = 4000
parallel = true

m<=Setting(:tpf_rstar, 2.0)
m<=Setting(:tpf_N_MH, N_MH)
m<=Setting(:tpf_c, 0.1)
m<=Setting(:tpf_accept_rate, 0.5)
m<=Setting(:tpf_target, 0.4)
m<=Setting(:tpf_n_particles, n_particles)
m<=Setting(:use_parallel_workers, parallel)
m<=Setting(:x_tolerance, zero(float(0)))
m<=Setting(:tpf_deterministic, true)

            df = readtable("$path/../../../test/reference/us.txt",header=false,separator=' ')
            data = convert(Matrix{Float64},df)
            data=data'

            # Read in matrices from Schorfheide Matlab code
            file = "$path/../../../test/reference/matlab_variable_for_testing.h5"
            A = h5read(file, "A")
            B = h5read(file, "B")
            H = h5read(file, "H")
            R = h5read(file, "R")
            S2 = h5read(file,"S2")
            Φ = h5read(file, "Phi")

rand_mat = randn(size(S2,1) ,1)
m<=Setting(:tpf_rand_mat, rand_mat)

	    #Write rand_mat to h5 for Matlab to read for comparison
            rand_mat = randn(size(S2,1),1)    	
	    h5open("$path/../../../test/reference/mutationRandomMatrix.h5","w") do file
                write(file,"rand_mat",rand_mat)
    	    end
            m<=Setting(:tpf_rand_mat,rand_mat)
            
            # Set variables within system
    	    transition_equation = Transition(Φ, R)
    	    measurement_equation = Measurement(B,squeeze(A,2),S2,H,rand_mat,R)
    	    system = System(transition_equation, measurement_equation)
            
            # Parameters given in Schorfheide's MATLAB code
            params = [2.09, 0.98, 2.25, 0.65, 0.34, 3.16, 0.51, 0.81, 0.98, 
                       0.93, 0.19, 0.65, 0.24,0.12,0.29,0.45]
            update!(m,params)


s0 = zeros(size(system[:TTT])[1])
P0 = nearestSPD(solve_discrete_lyapunov(Φ,R*S2*R'))

addprocs_sge(10,queue="background.q")
@everywhere using DSGE

neff, lik = tpf_matlab(m, data, system, s0, P0)

h5open("$path/../../../test/reference/output_likelihoods_ansch.h5","w") do f
     write(f,"julia_likelihoods",lik)
end

nothing