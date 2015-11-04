using Base.Test, Compat

function verbose_dict()
    @compat(Dict{Symbol,Int}(:none => 0, :low => 1, :high => 2))
end

function test_matrix_eq2{T<:AbstractFloat, U<:AbstractString}(expect::Array{T}, 
                                                              actual::Array{T}, 
                                                              expstr::U,
                                                              actstr::U,
                                                              ϵ_abs::Float64 = 1e-6,
                                                              ϵ_rel::Float64 = 1e-2)
    if length(expect) ≠ length(actual)
        error("lengths of ", expstr, " and ", actstr, " do not match: ",
              "\n  ", expstr, " (length $(length(expect))) = ", expect,
              "\n  ", actstr, " (length $(length(actual))) = ", actual)
    end

    # Absolute difference filter
    abs_diff   = actual .- expect .> ϵ_abs
    n_abs_diff = sum(abs_diff)

    # Relative difference filter
    rel_diff   = 100(actual .- expect) ./ expect .> ϵ_rel
    n_rel_diff = sum(rel_diff)

    # Element is only problematic if it fails *both* tests.
    mixed_diff   = abs_diff & rel_diff
    n_mixed_diff = sum(mixed_diff)

    if n_mixed_diff ≠ 0
        sdiff = string("|a - b| <= ", ϵ_abs, 
                   " or |a - b|/|b| <= ", ϵ_rel, "%,",
                   " ∀ a ∈ ", actstr, ",",
                   " ∀ b ∈ ",expstr)
        error("assertion failed:\n", 
              "    ", sdiff,
              "\n$(n_abs_diff) entries fail absolute equality filter",
              "\n$(n_rel_diff) entries fail relative equality filter",
              "\n$(n_mixed_diff) entries fail both equality filters")
    end

end

"""
    @test_matrix_approx_eq(a, b)

Test two matrices of floating point numbers `a` and `b` for approximate equality.
"""
macro test_matrix_approx_eq(a,b)
    :(test_matrix_eq2($(esc(a)),$(esc(b)),$(string(a)),$(string(b))))
end

"""
    @test_matrix_approx_eq_eps(a, b, ϵ_abs, ϵ_rel)

Test two matrices of floating point numbers `a` and `b` for equality taking in account a
margin of absolute tolerance given by `ϵ_abs` and a margin of relative tolerance given by
`ϵ_rel` (in comparison to `b`).
"""
macro test_matrix_approx_eq_eps(a,b,c,d)
    :(test_matrix_eq2($(esc(a)),$(esc(b)),$(string(a)),$(string(b)),$(esc(c)),$(esc(d))))
end

# minusnan(x, y) evaluates x-y in a way that treats NaN like Inf and sets Inf - Inf = 0
minusnan{T<:AbstractFloat}(x::T, y::T) =  minusnan(complex(x), complex(y))

function minusnan{T<:AbstractFloat}(x::Complex{T}, y::Complex{T})
    x = isnan(x) || isinf(x) ? Inf : x
    y = isnan(y) || isinf(y) ? Inf : y
    return isinf(x) && isinf(y) ? 0 : x - y
end

# percenterr(x,y) returns the absolute value of the percentage error of y wrt x
function percenterr{T<:AbstractFloat}(x::Complex{T}, y::Complex{T})
    return abs((x-y)/x)
end

function percenterr{T<:AbstractFloat}(x::T, y::T)
    return abs((x-y)/x)
end

# checksigns(x,y) checks to see if x and y have the same sign. returns 1 if they do, 0 if
# not, and -1 if one is zero but the other isnt
function checksigns(x::Number,y::Number)
    sx = sign(x)
    sy = sign(y)

    if sx == sy 
        return 1
    elseif isnan(x) && isnan(y)
        return 1
    elseif sx == 0 || sy == 0
        return -1
    else
        return 0
    end
end

# Compares matrices, reports absolute differences, returns true if all entries close enough
function test_matrix_eq{T<:AbstractFloat}(expect::Array{T}, 
                                          actual::Array{T}; 
                                          ϵ_abs::Float64 = 1e-3, 
                                          ϵ_rel::Float64 = 1e-2, 
                                          verbose::Bool = false)

    n_entries = length(expect)
    same_sign = map(checksigns, expect, actual)
    diff_sign_zeros = count(x-> x < 0, same_sign) # # that are different because 1 sign is 0
    n_diff_sign = count(x -> x != 1, same_sign)

    if verbose
        println("$n_diff_sign of $n_entries entries have opposite signs")
        println("$diff_sign_zeros of $n_entries entries have different signs, but one of the entries is 0")
    end

    abs_diff_eq = test_matrix_eq(complex(expect), complex(actual); 
    ϵ_abs=ϵ_abs, ϵ_rel=ϵ_rel, verbose=verbose)
    
    return abs_diff_eq 
end

# Complex-valued input matrices
function test_matrix_eq{T<:AbstractFloat}(expect::Array{Complex{T}},
                                          actual::Array{Complex{T}}; 
                                          ϵ_abs::Float64 = 1e-3, 
                                          ϵ_rel::Float64 = 1e-2, 
                                          verbose::Bool = false)

    # Matrices of different sizes return false
    if size(expect) != size(actual)
        if verbose
            println("Size expect $(size(expect)), actual $(size(actual))\n")
        end
        return false
    end

    n_dims = ndims(expect)
    n_entries = length(expect)

    # Count differences and find max
    abs_diff = abs(map(minusnan, expect, actual))
    n_neq = countnz(abs_diff)
    diff_inds = find(x -> x > ϵ_abs, abs_diff)    
    n_not_approx_eq = length(diff_inds)
    max_abs_diff = maximum(abs_diff)
    max_inds = ind2sub(size(abs_diff), indmax(abs_diff))

    # Count percent differences and find max
    pct_err = map(percenterr, expect, actual) #percenterr returns an absolute value
    n_neq_pct = countnz(pct_err)
    diff_inds_pct = find(x -> x > ϵ_rel, pct_err)
    n_not_approx_eq_pct = length(diff_inds_pct)
    max_pct_err = maximum(pct_err)
    max_pct_inds = ind2sub(size(pct_err), indmax(pct_err))

    # Find elements that are both absolutely and relatively different
    very_different = intersect(diff_inds, diff_inds_pct)
    
    # Print output
    if verbose
        println("$n_neq of $n_entries entries with abs diff > 0")
        println("$n_not_approx_eq of $n_entries entries with abs diff > $ϵ_abs")
        # println("$n_diff_sign of $n_entries entries have opposite signs")

        # If there are absolute differences
        if n_neq != 0
            println("Max abs diff of $max_abs_diff at entry $max_inds")
            if isa(expect, Matrix)
                println("The entries at $max_inds are $(expect[max_inds[1],max_inds[2]]) (expect) and $(actual[max_inds[1], max_inds[2]]) (actual)")
            end
            
        else 
            println("Max abs diff of 0")
        end

        # If there are percent differences
        if n_neq_pct != 0
            println("Max percent error of $max_pct_err at entry $max_pct_inds")
            if isa(expect, Matrix)
                println("The entries at $max_pct_inds are $(expect[max_pct_inds[1],max_pct_inds[2]]) (expect) and $(actual[max_pct_inds[1], max_pct_inds[2]]) (actual)\n")
            end
        else 
            println("Max pct diff of 0\n")
        end
    end

    # Warn if there are absolute diffs or relative diffs 
    if n_not_approx_eq !=0 
        warn("Absolute differences between entries exceed threshold of $ϵ_abs")
    end

    if n_not_approx_eq_pct !=0
        
        # Super-small numbers result in huge % differences between
        # numbers, so only warn if the numbers arent tiny
        for i = 1:length(diff_inds_pct)
            if real(expect[diff_inds_pct[i]]) > real(ϵ_rel)
                warn("Relative differences between entries exceed threshold of $ϵ_rel\%")
                break
            end
        end
    end
        
    # Fail if there are elements that have absolute differences > ϵ_abs and relative differences > ϵ_rel
    return length(very_different) == 0
end

# Tests whether the absolute values of 2 matrices are equal
function test_matrix_abs_eq{R<:AbstractFloat, S<:AbstractFloat, T<:AbstractFloat}(expect::
             Array{R}, actual::Array{S}; ϵ_abs::T = 1e-4, verbose::Bool = false)

    return test_matrix_eq(abs(expect), abs(actual), ϵ_abs=ϵ_abs, verbose=verbose)
end

function test_util()
    test_test_matrix_eq()
    test_readcsv_complex()
    test_test_eigs_eq()
    println("### All tests in util.jl tests passed\n")
end


function test_test_matrix_eq()
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
    @test test_matrix_eq(m1, m3; ϵ_abs=2e-4) # but true with larger ϵ_abs

    # Arguments of different float type returns true
    m2_float16 = convert(Matrix{Float16}, m2)
    ϵ_abs = convert(Float32, 0.0001)
    @test test_matrix_eq(m1, m2)

    # Complex-valued matrices
    @test test_matrix_eq(complex(m1), complex(m2))
    @test !test_matrix_eq(complex(m1), complex(m3))

    # 3D arrays
    m4 = zeros(2, 2, 2)
    m5 = ones(2, 2, 2)
    @test !test_matrix_eq(m4, m5)

    println("test_matrix_eq tests passed\n")
end

nothing
