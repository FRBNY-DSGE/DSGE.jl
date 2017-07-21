using ClusterManagers

parallel = true
T = 200

if parallel
    my_procs = addprocs_sge(10,queue="background.q")
    @everywhere using DSGE
end

time_blocks = zeros(T)

for t=1:T
    tic()
    parallel_breaking()
    print("Time $t ")
    time_blocks[t] = toc()
    print("===========================================")
end
rmprocs(my_procs)