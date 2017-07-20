
# Test to show buildup of parallel calls
function parallel_breaking()
    tic()
    out = @sync @parallel (hcat) for i=1:N
        move_stuff(A, B, C[:,i], D[:,i])
    end
    print("\nOne function call: ")
    toc()
end

function move_stuff(A::Array{Float64}, B::Array{Float64}, C::Array{Float64}, D::Array{Float64})
    
    
    return a, b, c
end 
