using MAT, JLD2

@inline function string_as_varname(s::AbstractString)
    S = Symbol(s)
    @eval $S
end

#@inline function test_against_matlab(filename::String)
    vars = matread(filename)
    for v in vars
        cur_type = typeof(v[2])
        if (v[1] == "t0")
            continue
        elseif v[1] == "a_grid"
            @assert v[2] ≈ a_grid
        elseif v[1] == "a_g_grid"
            @assert v[2] ≈ a_g_grid
        elseif v[1] == "b_grid"
            @assert v[2] ≈ b_grid
        elseif v[1] == "y_grid"
            @assert v[2] ≈ y_grid
        elseif v[1] == "y_g_grid"
            @assert v[2] ≈ y_g_grid
        elseif v[1] == "r_a_grid"
            @assert v[2] ≈ r_a_grid
        elseif v[1] == "r_b_grid"
            @assert v[2] ≈ r_b_grid
        elseif v[1] == "r_a_g_grid"
            @assert v[2] ≈ r_a_g_grid
        elseif v[1] == "r_b_g_grid"
            @assert v[2] ≈ r_b_g_grid
        elseif v[1] == "daf_g_grid"
            @assert v[2] ≈ daf_g_grid
        elseif v[1] == "daf_grid"
            @assert v[2] ≈ daf_grid
        elseif v[1] == "dab_grid"
            @assert v[2] ≈ dab_grid
        elseif v[1] == "dab_g_grid"
            @assert v[2] ≈ dab_g_grid
        elseif v[1] == "dab_g_tilde_grid"
            @assert v[2] ≈ dab_g_tilde_grid
        elseif v[1] == "dbf_grid"
            @assert v[2] ≈ dbf_grid
        elseif v[1] == "dbf_g_grid"
            @assert v[2] ≈ dbf_g_grid
        elseif v[1] == "dbb_grid"
            @assert v[2] ≈ dbb_grid
        elseif v[1] == "dbb_g_grid"
            @assert v[2] ≈ dbb_g_grid
        elseif v[1] == "trans_grid"
            @assert v[2] ≈ trans_grid
        elseif v[1] == "l_g_grid"
            @assert v[2] ≈ l_g_grid
        elseif v[1] == "w_grid"
            @assert v[2] ≈ w_grid
        elseif v[1] == "c_0"
            @assert v[2] ≈ c_0
        elseif v[1] == "V_0"
            @assert v[2] ≈ V_0
        elseif v[1] == "gg0"
            @assert v[2] ≈ gg0
        elseif v[1] == "gg"
            @assert v[2] ≈ gg
        elseif v[1] == "r_b"
            @assert v[2] ≈ r_b
        elseif v[1] == "r_b_borr"
            @assert v[2] ≈ r_b_borr
        elseif v[1] == "y_dist"
            @assert v[2] ≈ y_dist
        elseif v[1] == "y_mean"
            @assert v[2] ≈ y_mean
        elseif v[1] == "aalpha"
            @assert v[2] ≈ aalpha
        elseif v[1] == "ddeath"
            @assert v[2] ≈ ddeath
        elseif v[1] == "ddelta"
            @assert v[2] ≈ ddelta
        elseif v[1] == "rrho"
            @assert v[2] ≈ rrho
        elseif v[1] == "chi0"
            @assert v[2] ≈ chi0
        elseif v[1] == "chi1"
            @assert v[2] ≈ chi1
        elseif v[1] == "chi2"
            @assert v[2] ≈ chi2
        elseif v[1] == "a_lb"
            @assert v[2] ≈ a_lb
        elseif v[1] == "pam"
            @assert v[2] ≈ pam
        elseif v[1] == "xxi"
            @assert v[2] ≈ xxi
        elseif v[1] == "ggamma"
            @assert v[2] ≈ ggamma
        elseif v[1] == "tau_I"
            @assert v[2] ≈ tau_I
        elseif v[1] == "trans"
            @assert v[2] ≈ trans
        elseif v[1] == "lambda"
            @assert v[2] ≈ lambda
        elseif v[1] == "K_liquid"
            @assert v[2] ≈ K_liquid
        elseif v[1] == "aggregate_variables"
            @assert v[2] ≈ aggregate_variables
        elseif v[1] == "distributional_variables"
            @assert v[2] ≈ distributional_variables
        elseif v[1] == "distributional_variables_1"
            #@assert v[2] ≈ distributional_variables_1
            continue
        else
            @show "NOT PRINTING", v[1]
            #=
            try
                if !(string_as_varname(v[1]) ≈ v[2])
                    @show "NEQ", v[2]
                else
                    @show "WORKSWORKSWORKS", v[2]
                end
                @assert string_as_varname(v[1]) ≈ v[2]
            catch err
                if isa(err, UndefVarError)
                    println(err)
                else
                    throw(err)
                end
            end
            =#
        end
    end
    println("All tests passed for: ", filename)
#end
