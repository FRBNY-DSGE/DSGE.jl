#= The following script creates 30 workers. The function problem_short sets up inputs then for T time steps (default 50) runs the mutation function 10 times per time step. It keeps track of runtimes and prints a plot of the increase in runtime for each time step. It calls mutation_problem.jl instead of mutation.jl which are identical however mutation_problem.jl has a flag for whether to create two MvNormal objects
=#

#= Two things that we've done that seem to mitigate the problem:
1. Not passing in the model object (relatively large) into mutation--> this would 
seem to suggest the problem is with allocating work
2. Not creating MvNormal Distribution objects in mutation--> this occurs within
 each mutation which would seem to suggest the problem is with gunk accumulating on each worker.

Also, the problem is worse with 30 workers than 10 workers so to see an effect, should run with
30 workers

Also, look at distVsNodist.pdf in the src/estimate folder
=#

using ClusterManagers, HDF5, Plots, JLD
using QuantEcon: solve_discrete_lyapunov

addprocs_sge(30, queue="background.q")

@everywhere using DSGE

function the_problem(parallel=true::Bool,T=50::Int64, distCall=true)
    ## Setup
    m,system,TTT,sqrtS2,s0,P0,s_lag_tempered,ε,yt,nonmissing, N_MH,c, n_particles,deterministic, μ, cov_s,s_t_nontempered = setup()
   
    ## Initialize vector for keeping track of runtimes for each time step
    times = zeros(T)
    
    # Loop over number of desired time steps
    for t = 1:T
        @show t
        #Begin timer for time step
        tic()        
        #Run mutation 10 times per time step
        for i=1:10
            acpt_vec=zeros(n_particles)
            print("Mutation ")
            #Begin timer for each mutation step
            tic()
            if parallel
                print("(in parallel) ")
                #out = pmap(i->mutation_problem(c, N_MH,deterministic,yt,s_lag_tempered[:,i],ε[:,i],cov_s,nonmissing,distCall), 1:n_particles)
                out = @sync @parallel (hcat) for i=1:n_particles
                    mutation_problem(c,N_MH,deterministic,system,yt,s_lag_tempered[:,i],ε[:,i],cov_s,nonmissing,distCall)
                end
            else
                print("(not parallel)")
                out = [mutation_problem(c,N_MH,deterministic,system,yt,s_lag_tempered[:,i],ε[:,i],cov_s,nonmissing,distCall) for i=1:n_particles]
            end               
            toc()
            # Disentangle three outputs of mutation and enter them into appropriate arrays
            for i = 1:n_particles
                s_t_nontempered[:,i] = out[i][1]
                ε[:,i] = out[i][2]
                acpt_vec[i]=out[i][3]
            end
        end
        print("Completion of one period ")
        gc()
        times[t] = toc()
    end
    h5open("../../test/reference/the_problem_times_distCall=$distCall.h5", "w") do file
        write(file, "times", times)
    end    
    plotly()
    plot(times)
    gui()
end


"""
```
get_chol(mat::Aray)
```
Calculate and return the Cholesky of a matrix.
"""
function get_chol(mat::Array{Float64,2})
    return Matrix(chol(nearestSPD(mat)))
end

function setup()
    m = SmetsWouters("ss1", testing=true)

    path = dirname(@__FILE__)
    filesw= "/data/dsge_data_dir/dsgejl/realtime/input_data/data"
    data = readcsv("$filesw/realtime_spec=smets_wouters_hp=true_vint=110110.csv",header=true)
    data = convert(Array{Float64,2},data[1][:,2:end])
    data=data'
 
    N_MH= 10
    c = 0.1
    n_particles = 4000

    system = load( "$path/../../test/reference/system.jld", "system")
    RRR = system.transition.RRR
    TTT = system.transition.TTT
    S2 = system.measurement.QQ
    sqrtS2 = RRR*get_chol(S2)'
   
    s0 = zeros(size(TTT)[1])
    P0 = nearestSPD(solve_discrete_lyapunov(TTT,RRR*S2*RRR'))

    n_errors = size(S2,1)
    n_states = size(system.measurement.ZZ,2)

    s_lag_tempered_rand_mat = randn(n_states,n_particles)
    ε = randn(n_errors, n_particles)
   
    s_lag_tempered = repmat(s0,1,n_particles) + get_chol(P0)'*s_lag_tempered_rand_mat
    
    yt = data[:,25]
    nonmissing =!isnan(yt)
    deterministic=false

    μ = mean(ε,2)
    cov_s = (1/n_particles)*(ε-repmat(μ,1,n_particles))*(ε - repmat(μ,1,n_particles))'
    if !isposdef(cov_s)
        cov_s = diagm(diag(cov_s))
    end
    
    s_t_nontempered = TTT*s_lag_tempered + sqrtS2*ε
      
    return m,system,TTT,sqrtS2,s0,P0,s_lag_tempered,ε,yt,nonmissing, N_MH,c, n_particles, deterministic, μ, cov_s, s_t_nontempered
end