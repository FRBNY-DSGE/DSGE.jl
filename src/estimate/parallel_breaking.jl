
# Test to show buildup of parallel calls
function parallel_breaking()
    N = 4000
    K = 7
    T = 100
    m=Array{Array{Float64}}(4)
    time_blocks = zeros(T)
    A = randn(N,K)
    B = randn(K,N)
    C = randn(K,N)
    D = randn(N,N)
    j = 1

    m[1]=A
    m[2]=B
    m[3]=C
    m[4]=D

    for t=1:1
        tic()
        while j<5
            tic()
            out = @allocated @sync @parallel (hcat) for i=1:N
                move_stuff(m, A, B, C[:,i], D[:,i])
            end
            print("\nOne function call ")
            toc()
            j+= 1
        end
        print("Time $t ")
        time_blocks[t] = toc()
        print("==========================")
    end
    return time_blocks
end

function move_stuff(m::Array, A::Array{Float64}, B::Array{Float64}, C::Array{Float64}, D::Array{Float64})
    A = m[1]
    B = m[2]
    C = m[3]
    D = m[4]

    M = A*C + D
    a = b = c = 0
    for i=1:3
        B = B'
        M1 = A*C + D
        M2 = A*C + D
    end
    gc()
    return a, b, c
    
end 
