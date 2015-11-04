using Base.Test, Compat

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

nothing
