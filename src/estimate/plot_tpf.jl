function plot_tpf()
#using ClusterManagers, HDF5, DataFrames, Plots, StateSpaceRoutines
#using QuantEcon: solve_discrete_lyapunov


  #  custom_settings= Dict{Symbol,Setting}(
#        :date_forecast_start => Setting(:date_forecast_start, quartertodate("2011-Q2")))
    m = SmetsWouters("ss1", testing=true)
#    m <= Setting(:date_conditional_end,quartertodate("2011-Q1"))
 #   m <= Setting(:date_forecast_start, quartertodate("2011-Q2"))
    
    path = dirname(@__FILE__)
    
    filesw = "/data/dsge_data_dir/dsgejl/realtime/input_data/data"
    data = readcsv("$filesw/realtime_spec=smets_wouters_hp=true_vint=110110.csv",header=true)
    data=convert(Array{Float64,2},data[1][:,2:end])
    data = data'
    
    params = h5read("$filesw/../../output_data/smets_wouters/ss0/estimate/raw/paramsmode_vint=110110.h5","params")
    push!(params, m[:e_y].value, m[:e_L].value, m[:e_w].value, m[:e_π].value, m[:e_R].value, m[:e_c].value, m[:e_i].value)
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

#For comparing change variance to not changing variance    
    Neff, lik, times = tpf(m,data,system,s0,P0,false)
    kalman_out = filter(m,data,s0,P0,allout=true)
    NeffVar, likVar, timesVar = tpf(m,data,system,s0,P0,true)

    h5open("$path/../../test/reference/varLik.h5","w") do file
        write(file, "lik", lik)
        write(file, "likVar",likVar)
        write(file, "kalman_lik",kalman_out[:L_vec])
    end

    
##For getting distribution of errors (relative to Kalman filter)
    # kalman_out = filter(m,data,s0,P0,allout=true)
    # for i=1:2
    #     @show i
    #     Neff, lik = tpf(m,data,system,s0,P0,true)
    #     error[i] = kalman_out[:L]-sum(lik)
    #     h5open("$path/../../test/reference/error.h5","w") do file
    #         write(file, "error", error)
    #     end
    # end

    # plotly()
    # plot(kalman_out, color=:blue, linewidth=3, label="Kalman")
    # plot!(lik, color=:green, linewidth=2, label="TPF")
    # plot!(legend=:bottomright, xlabel="Time", ylabel="Log likelihood")
    # gui()
    
end
    #kalman filter
    #=kalman_out = filter(m,data,s0,P0,allout=true)
    kalman_out[:L_vec]

    neff, lik = tpf(m,data,system,s0,P0)
    h5open("$path/../../test/reference/lik_for_plot_4000_10mh.h5","w") do file
        write(file, "kalman_lik",kalman_out[:L_vec])
        write(file, "tpf_lik", lik)
    end
    @show mean(abs(lik-kalman_lik))
=#
#=
    #for getting distribution of errors
    kalman_out = filter(m,data,s0,P0,allout=true)
    error = zeros(50)

    out = pmap(i-> tpf(m,data,system,s0,P0), 1:100)
    for i = 1:100
        error[i]=kalman_out[:L]-sum(out[2])
    end

=#



###########
#=    
    srand(1234)
    df = readtable("$path/../../test/reference/us.txt",header=false,separator=' ')
    data=convert(Matrix{Float64},df)
    data=data'
    file = "$path/../../test/reference/matlab_variable_for_testing.h5"
    A = h5read(file,"A")
    B = h5read(file, "B")
    H = h5read(file, "H")
    R = h5read(file, "R")
    S2 = h5read(file, "S2")
    Φ = h5read(file, "Phi")
    m<=Setting(:DD,A[:,1])
    m<=Setting(:ZZ,B)
    m<=Setting(:RRR,R)
    m<=Setting(:TTT,Φ)
    m<=Setting(:EE,H)
    m<=Setting(:tpf_S2,S2)
    C= zeros(8,1)
    m<=Setting(:CCC,C)
    params=[2.09,0.98,2.25,0.65,0.34,3.16,0.51,0.81,0.98,0.93,0.19,0.65,0.24,0.12,0.29,0.45]
    update!(m,params)
    s0=zeros(8)
    P0 = nearestSPD(solve_discrete_lyapunov(Φ, R*S2*R'))
  
  h5open("$path/../../test/reference/output_likelihoods_500_testing.h5","w") do file
        write(file, "julia_likelihoods",lik)
    end
   
 
    #4000 particles. testing
    n_particles = 4000
    m<=Setting(:tpf_n_particles,n_particles)
    neff,lik = tpf(m,data,s0,P0)
    h5open("$path/../../test/reference/output_likelihoods_4000_testing.h5","w") do file
        write(file,"julia_likelihoods",lik)
    end
    
    #4000 particles. no testing]
    df = readtable("$path/../../test/reference/us.txt", header=false, separator=' ')
    data = convert(Matrix{Float64},df)
    data=data'
    
    file = "$path/../../test/reference/matlab_variable_for_testing.h5"
    A = h5open(file, "r") do file
        read(file,"A")
    end
    B = h5open(file,"r") do file
        read(file,"B")
    end
    H = h5open(file,"r") do file
        read(file,"H")
    end
    R = h5open(file,"r") do file
        read(file,"R")
    end
    S2 = h5open(file,"r") do file
        read(file,"S2")
    end
    Φ = h5open(file,"r") do file
        read(file, "Phi")
    end
    m<=Setting(:DD,A[:,1])
    m<=Setting(:ZZ,B)
    m<=Setting(:RRR,R)
    m<=Setting(:TTT,Φ)
    m<=Setting(:EE,H)
    m<=Setting(:tpf_S2,S2)
    C = zeros(8,1)
    m<=Setting(:CCC,C)
    params=[2.09,0.98,2.25,0.65,0.34,3.16,0.51,0.81,0.98,0.93,0.19,0.65,0.24,0.12,0.29,0.45]
    update!(m,params)
    
    neff,lik = tpf(m,data,s0,P0,testing=false)
    h5open("$path/../../test/reference/output_likelihoods_4000_notesting.h5","w") do file
        write(file,"julia_likelihoods",lik)
    end
    
    #Plot run-time vs. number of particles
    sys = compute_system(m)
    R = sys.transition.RRR
    S2 = sys.measurement.QQ
    Φ = sys.transition.TTT
    file = "$path/../../test/reference/optimize.h5"
    x0 = h5read(file,"params")
    data = h5read(file,"data")
    miniimzer = h5read(file,"data")
    minimizer = h5read(file, "minimizer")
    update!(m,x0)
    x0 = Float64[p.value for p in m.parameters]
    out,H = optimize!(m,data;iterations=500)
    params = out.minimizer
    @show params
    update!(m,params)
    s0 = zeros(8)
    P0 = nearestSPD(solve_discrete_lyapunov(Φ, R*S2*R'))
                    
    times = []
    ns = [100,300,500,700,900,1100,1300,1500,1700,1900,2100,4000]
    
    times = pmap(n->time(m,data,s0,P0,false,n),ns)
    @show times
    plotly()
    plot(ns,times)
    gui()
end
=#