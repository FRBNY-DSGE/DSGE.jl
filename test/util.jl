# Utilities for testing DSGE.jl



# TODO: decide what a sensible default ε value is for our situation
# Compares matrices, reports absolute differences, returns true if all entries close enough
function test_matrix_eq{R<:FloatingPoint, S<:FloatingPoint, T<:FloatingPoint}(expected::Array{R, 2},
             actual::Array{S, 2}, ε::T = 1e-4)
    # Matrices of different sizes return false
    if size(expected) != size(actual)
        return false
    end

    # Size variables and counters
    (rows, cols) = size(expected)
    n_entries = rows * cols
    
    n_neq = 0                   # number of entries with abs diff > 0
    n_not_approx_eq = 0         # number of entries with abs diff > ε
    max_abs_diff = 0.0          # maximum abs diff
    max_inds = (0, 0)           # indices of maximum abs diff

    # Count differences and find max
    for j = 1:cols
        for i = 1:rows
            abs_diff = abs(expected[i, j] - actual[i, j])
            if abs_diff > 0
                n_neq += 1
                if abs_diff > ε
                    n_not_approx_eq += 1
                end
                if abs_diff > max_abs_diff
                    max_abs_diff = abs_diff
                    max_inds = (i, j)
                end
            end
        end
    end

    # Print output
    println("$n_neq of $n_entries entries with abs diff > 0")
    println("$n_not_approx_eq of $n_entries entries with abs diff > $ε")
    if max_inds != (0, 0)
        (i, j) = max_inds
        println("Max abs diff of $max_abs_diff at entry ($i, $j)\n")
    else
        println("Max abs diff of 0\n")
    end

    # Return true if all entries within ε
    return n_not_approx_eq == 0
end



function test_util()
    # Matrices of different sizes returns false
    m0 = zeros(2, 3)
    m1 = zeros(2, 2)
    @test !test_matrix_eq(m0, m1)

    # Returns true
    m2 = [0.0001 0.0; 0.0 -0.0001]
    @test test_matrix_eq(m1, m2)

    # Returns false
    m3 = [0.0001 0.0; 0.0 -0.0002]
    @test !test_matrix_eq(m1, m3)
    @test test_matrix_eq(m1, m3, 2e-4) # but true with larger ε

    # Arguments of different float type returns true
    m2_float16 = convert(Array{Float16, 2}, m2)
    ε = convert(Float32, 0.0001)
    @test test_matrix_eq(m1, m2)

    println("util.jl tests passed")
end
