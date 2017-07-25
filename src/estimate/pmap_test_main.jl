function time1(parallel::Bool, pmap_flag::Bool)
    a = zeros(500)
    b = zeros(500)
    c = zeros(500)
    if pmap_flag & !parallel
        tic()
        arr = pmap(x->pmap_test1(x),1:500)
        for i = 1:500
            a[i]=arr[i][1]
            b[i]=arr[i][2]
            c[i]=arr[i][3]
        end
        toc()
    elseif parallel & !pmap_flag
        tic()
        arr = @sync @parallel (hcat) for x=1:500
            pmap_test1(x)
        end
        for i = 1:500
            a[i]=arr[i][1]
            b[i]=arr[i][2]
            c[i]=arr[i][3]
        end
        toc()
    else
        tic()
        arr = [pmap_test1(x) for x = 1:500]
        for i = 1:500
            a[i]=arr[i][1]
            b[i]=arr[i][2]
            c[i]=arr[i][3]
        end
        toc()
    end
end

function time2(parallel::Bool, pmap_flag::Bool)
    a = zeros(500)
    b = zeros(500)
    c = zeros(500)
    A = randn(500,8)
    B = randn(8,500)
    C = randn(1000,500)
    D = randn(500,700)
    if pmap_flag & !parallel
        tic()
        arr = pmap(x->pmap_test2(A,B,C,D),1:500)
        for i = 1:500
            a[i]=arr[i][1]
            b[i]=arr[i][2]
            c[i]=arr[i][3]
        end
        toc()
    elseif parallel & !pmap_flag
        tic()
        arr = @sync @parallel (hcat) for x=1:500
            pmap_test2(A,B,C,D)
        end
        for i = 1:500
            a[i]=arr[i][1]
            b[i]=arr[i][2]
            c[i]=arr[i][3]
        end
        toc()
    else
        tic()
        arr = [pmap_test2(A,B,C,D) for x = 1:500]
        for i = 1:500
            a[i]=arr[i][1]
            b[i]=arr[i][2]
            c[i]=arr[i][3]
        end
        toc()
    end
end


function time3(parallel::Bool, pmap_flag::Bool)
    a = zeros(500)
    b = zeros(500)
    c = zeros(500)
    A = randn(500,8)
    B = randn(8,500)
    C = randn(1000,500)
    D = randn(500,700)
    if pmap_flag & !parallel
        tic()
        arr = pmap(x->pmap_test3(A,B,C,D),1:500)
        for i = 1:500
            a[i]=arr[i][1]
            b[i]=arr[i][2]
            c[i]=arr[i][3]
        end
        toc()
    elseif parallel & !pmap_flag
        tic()
        arr = @sync @parallel (hcat) for x=1:500
            pmap_test3(A,B,C,D)
        end
        for i = 1:500
            a[i]=arr[i][1]
            b[i]=arr[i][2]
            c[i]=arr[i][3]
        end
        toc()
    else
        tic()
        arr = [pmap_test3(A,B,C,D) for x = 1:500]
        for i = 1:500
            a[i]=arr[i][1]
            b[i]=arr[i][2]
            c[i]=arr[i][3]
        end
        toc()
    end
end


function time4(parallel::Bool, pmap_flag::Bool)
    a = [zeros(3) for i=1:4000]
    b = [zeros(3) for i=1:4000]
    c = zeros(4000)
    m = SmetsWouters("ss1",testing=true)
    yt = randn(7)
    # D = randn(7)
    # Z = randn(7,54)
    # E = diagm(rand(7))
    # R = randn(54,7)
    # T = randn(54,54)
    # Q = diagm(rand(7))
    # sqrtS2 = R*Matrix(chol(nearestSPD(Q)))'
    system = compute_system(m)
    nonmissing = [true, true, true, true, true, true, true] 
   
    @show nonmissing
    s_init = randn(54)
    ε_init = randn(7)
    cov_s = diagm(rand(7))

    if pmap_flag & !parallel
        tic()
        arr = pmap(x->pmap_test4(yt,system,s_init,ε_init, cov_s,nonmissing),1:4000)
        for i = 1:4000
            a[i]=arr[i][1]
            b[i]=arr[i][2]
            c[i]=arr[i][3]
        end
        toc()
    elseif parallel & !pmap_flag
        tic()
        arr = @sync @parallel (hcat) for x=1:4000
            pmap_test4(yt,system,s_init,ε_init, cov_s,nonmissing)
        end
        
        for i = 1:4000
            a[i]=arr[i][1]
            b[i]=arr[i][2]
            c[i]=arr[i][3]
        end
        toc()
    else
        tic()
        arr = [pmap_test4(yt,system,s_init,ε_init, cov_s,nonmissing) for x = 1:4000]
        for i = 1:4000
            a[i]=arr[i][1]
            b[i]=arr[i][2]
            c[i]=arr[i][3]
        end
        toc()
    end
end

function time5()
    m = SmetsWouters("ss1",testing=true)
    A = randn(100,100)
    B = randn(100,100) 
    x = randn(4000)
    tic()
    arr = @sync @parallel (hcat) for i=1:4000
        pmap_test5(m,A,B,x[i])
    end
    return toc()
end


using ClusterManagers, Plots, HDF5

addprocs_sge(30,queue="background.q")
@everywhere include("pmap_test.jl")
 
# println("@parallel on")
# time2(true,false)
# println("pmap on")
# time2(false,true)
# rmprocs()
# println("no parallel")
# time2(false,false)

@show time5()
 
times = zeros(3000)
for i = 1:3000
    @show i
    println("@parallel on")
    times[i] = time5()
   # println("pmap on")
    #time4(false,true)
end
path = dirname(@__FILE__)
h5open("$path/../../test/reference/timeOverTime.h5","w") do file
    write(file,"times", times)
end

plotly()
plot(times)
gui()

rmprocs()
aaa
println("Non-parallel")
time4(false, false)