using ClusterManagers

parallel = true

if parallel
    my_procs = addprocs_sge(10,queue="background.q")
    @everywhere using DSGE
    @everywhere using Base
end
 
times = @allocated parallel_breaking()
open("parallel_time_increases.txt","w") do f
    print(f, times)
end

rmprocs(my_procs)