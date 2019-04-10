using Statistics, SparseArrays
"""
``
eqcond(m::TwoAssetHANK)
```

Expresses the equilibrium conditions in canonical form using Γ0, Γ1, C, Ψ, and Π matrices.

### Outputs

* `Γ0` (`n_states` x `n_states`) holds coefficients of current time states.
* `Γ1` (`n_states` x `n_states`) holds coefficients of lagged states.
* `C`  (`n_states` x `1`) is a vector of constants
* `Ψ`  (`n_states` x `n_shocks_exogenous`) holds coefficients of iid shocks.
* `Π`  (`n_states` x `n_states_expectational`) holds coefficients of expectational states.
"""
function eqcond_lite(m::TwoAssetHANK)
    # Read in parameters
    aalpha = m[:aalpha].value::Float64
    ddelta = m[:ddelta].value::Float64
    ddeath = m[:ddeath].value::Float64
    rrho   = m[:rrho].value::Float64
    chi0   = m[:chi0].value::Float64
    chi1   = m[:chi1].value::Float64
    chi2   = m[:chi2].value::Float64
    a_lb   = m[:a_lb].value::Float64
    kappa  = m[:kappa].value::Float64
    pam    = m[:pam].value::Float64
    xxi    = m[:xxi].value::Float64
    ggamma = m[:ggamma].value::Float64
    tau_I  = m[:tau_I].value::Float64
    trans  = m[:trans].value::Float64
    n_SS   = m[:n_SS].value::Float64

    nnu_aggZ    = m[:nnu_aggZ].value::Float64
    ssigma_aggZ = m[:ssigma_aggZ].value::Float64

    # Set liquid rates
    r_b_SS       = m[:r_b_SS].value::Float64
    r_b_borr_SS  = m[:r_b_borr_SS].value::Float64
    borrwedge_SS = m[:borrwedge_SS].value::Float64

    #lambda                     = get_setting(m, :lambda)::Matrix{Float64}

    K_liquid                   = get_setting(m, :K_liquid) ? 1 : 0::Int64
    aggregate_variables        = get_setting(m, :aggregate_variables)::Int64
    distributional_variables   = get_setting(m, :distributional_variables)::Int64
    distributional_variables_1 = get_setting(m, :distributional_variables_1)::Int64
    permanent                  = get_setting(m, :permanent) ? 1 : 0::Int64

    I      = get_setting(m, :I)::Int64
    J      = get_setting(m, :J)::Int64
    I_g    = get_setting(m, :I_g)::Int64
    J_g    = get_setting(m, :J_g)::Int64
    N      = get_setting(m, :N)::Int64

    #a      = get_setting(m, :a)::Vector{Float64}
    #b      = get_setting(m, :b)::Vector{Float64}
    #a_g    = get_setting(m, :a_g)::Vector{Float64}
    #b_g    = get_setting(m, :b_g)::Vector{Float64}

    agrid_new = get_setting(m, :agrid_new)::Int64
    bgrid_new = get_setting(m, :bgrid_new)::Int64
    ygrid_new = get_setting(m, :ygrid_new)::Int64

    amin = get_setting(m, :amin)
    bmin = get_setting(m, :bmin)
    amax = get_setting(m, :amax)
    bmax = get_setting(m, :bmax)

    #y      = vec(get_setting(m, :y))::Vector{Complex{Float64}}
    #y_dist = get_setting(m, :y_dist)::Vector{Complex{Float64}}
    #y_mean = get_setting(m, :y_mean)::Complex{Float64}
    KL     = get_setting(m, :KL_0)::Float64
    r_b_fix= get_setting(m, :r_b_fix) ? 1 : 0::Int64

    n_v = get_setting(m, :n_v)::Int64
    n_g = get_setting(m, :n_g)::Int64
    n_p = get_setting(m, :n_p)::Int64
    n_Z = get_setting(m, :n_Z)::Int64
    nVars    = get_setting(m, :nVars)::Int64
    nEErrors = get_setting(m, :nEErrors)::Int64

    # Set liquid rates
    r_b_borr = r_b_borr_SS

    vars_SS     = vec(m[:vars_SS].value)::Vector{Float64}
    V_SS = vars_SS[1:n_v]
    g_SS = vars_SS[n_v + 1 : n_v + n_g]
    K_SS = vars_SS[n_v + n_g + 1]
    r_b_SS = vars_SS[n_v + n_g + 2]
    if aggregate_variables == 1
        aggY_SS = vars_SS[n_v+n_g+3] # aggregate output
        aggC_SS = vars_SS[n_v+n_g+4] # aggregate consumption
    elseif distributional_variables == 1
        C_Var_SS = vars_SS[n_v+n_g+3] # consumption inequality
        earn_Var_SS = vars_SS[n_v+n_g+4] # earnings inequality
    elseif distributional_variables_1 == 1
        C_WHTM_SS = vars_SS[n_v+n_g+3] # consumption of wealthy hand-to-mouth
        C_PHTM_SS = vars_SS[n_v+n_g+4] # consumption of poor hand-to-mouth
    end
    aggZ_SS = vars_SS[n_v+n_g+n_p+1] # aggregate Z

####

    a, a_g, a_g_0pos          = create_a_grid(agrid_new, J, J_g, amin, amax)
    _, _, _, b, b_g, b_g_0pos = create_b_grid(bgrid_new, I, I_g)
    lambda, y, y_mean, y_dist, _ = create_y_grid(N, ygrid_new)

    indx = 1
    V_SS = reshape(V_SS, I*J, N)
    g_SS_ind = (indx != N) ? g_SS[(indx - 1)*I_g*J_g + 1 : indx*I_g*J_g] :
        vcat(g_SS[(indx - 1)*I_g*J_g + 1 : end], 0.0)

    n_v = Int(n_v/N)
    n_g = Int((n_g + 1)/N)
    nVars = n_v + n_g + n_p + 1

    @inline function get_residuals_lite(vars::Vector{T}) where {T<:Real}
        # ------- Unpack variables -------
        V   = reshape(vars[1:n_v] .+ V_SS[indx], I, J)  # value function
        g   = vars[n_v + 1 : n_v + n_g] .+ g_SS_ind     # distribution
        K   = vars[n_v + n_g + 1] + K_SS                # aggregate capital
        r_b = vars[n_v + n_g + 2] + r_b_SS

        if aggregate_variables == 1
            aggY     = vars[n_v+n_g+3] + aggY_SS     # aggregate output
            aggC     = vars[n_v+n_g+4] + aggC_SS     # aggregate consumption
        elseif distributional_variables == 1
            C_Var    = vars[n_v+n_g+3] + C_Var_SS    # consumption inequality
            earn_Var = vars[n_v+n_g+4] + earn_Var_SS # earnings inequality
        elseif distributional_variables_1 == 1
            C_WHTM  = vars[n_v+n_g+3] + C_WHTM_SS    # consumption of wealthy hand-to-mouth
            C_PHTM  = vars[n_v+n_g+4] + C_PHTM_SS    # consumption of poor hand-to-mouth
        end
        aggZ       = vars[n_v+n_g+n_p+1] + aggZ_SS   # aggregate Z
        V_Dot      = vars[nVars + 1 : nVars + n_v]
        g_Dot      = vars[nVars + n_v + 1:nVars + n_v + n_g]
        aggZ_Dot   = vars[nVars + n_v + n_g + n_p + 1]
        VEErrors   = vars[2*nVars + 1 : 2 * nVars + n_v]
        aggZ_Shock = vars[2*nVars + nEErrors + 1]
        # ------- Unpack variables -------

        dab_tilde_grid, dab_g_tilde_grid, dab_g_tilde_mat, dab_g_tilde = set_grids_lite2(a,
                                                                           b, a_g, b_g)
        a_gg = repeat(a_g, inner=I_g)
        b_gg = repeat(b_g, outer=J_g)

        # Construct problem functions
        util, deposit, cost = construct_problem_functions(ggamma, chi0, chi1, chi2, a_lb)

        r_b_vec, r_b_g_vec, daf_vec, daf_g_vec, dab_vec, dab_g_vec, dab_tilde, dab_g_tilde, dbf_vec, dbf_g_vec, dbb_vec, dbb_g_vec, dab, dab_tilde_mat, dab_g, dab_g_tilde_mat = set_vectors_lite(a, b, a_g, b_g, r_b, r_b_borr)

        dab_aux   = reshape(dab, I*J, 1)
        dab_g_aux = reshape(dab_g, I_g*J_g, 1)

        loc = findall(b .== 0)
        dab_g_tilde_mat_inv = spdiagm(0 => vec(1.0 ./ dab_g_tilde))
        dab_g_small = reshape(dab_g, I_g * J_g, 1)

        dab_g_small           = dab_g_small ./ dab_g_small[loc] * ddeath
        dab_g_small[loc]     .= 0.0
        death_process         = -ddeath * my_speye(I_g * J_g)
        death_process[loc,:]  = vec(dab_g_small)
        #death_process         = kron(my_speye(N), death_process)

        daba_g_aux = dab_g_aux .* a_gg
        dabb_g_aux = dab_g_aux .* b_gg

        if indx == N
            g_end = (1 - sum(g .* dab_g_aux[1:end-1])) / dab_g[I_g, J_g]
            @show g_end
            gg    = vcat(g, g_end)
        end
        # Prices
        w   = (1 - aalpha) * (K ^ aalpha) * n_SS ^ (-aalpha) *
            ((permanent == 0) ? exp(aggZ) ^ (1-aalpha) : 1.)
        r_a = aalpha * (K ^ (aalpha - 1)) * (((permanent == 0) ? exp(aggZ) : 1.) *
            n_SS) ^ (1 - aalpha) - ddelta

        # Auxiliary variables
        r_b_borr = r_b .+ borrwedge_SS

        # SET GRIDS
        r_b_vec = r_b .* (b .>= 0) + r_b_borr .* (b .< 0)

        # Other necessary objects
        y_shock      = y .* exp.(kappa * aggZ * (y .- y_mean) ./ std(y))
        y_shock_mean = dot(y_shock, y_dist)
        y_shock      = real(y_shock ./ y_shock_mean .* y_mean)

        # ripped out
        println("Timing: solve_hjb()")
        @time c, s, d  = solve_hjb_lite(V, I_g, J_g, a_lb, ggamma, permanent,
                                   ddeath, pam, aggZ, xxi, tau_I, w, trans,
                                   r_b_vec, y_shock[indx], a, b, cost, util, deposit)

        interp_decision = interpTwoD(b_g, a_g, b, a)#kron(my_speye(N), interpTwoD(b_g, a_g, b, a))
        d_g = reshape(interp_decision * vec(d), I_g, J_g)
        s_g = reshape(interp_decision * vec(s), I_g, J_g)
        c_g = reshape(interp_decision * vec(c), I_g, J_g)

        # Derive transition matrices
        println("Timing: transition_deriva()")
        @time A, AT = transition_deriva_lite(permanent==1, ddeath, pam, xxi, w, a_lb, aggZ,
                                        d, d_g, s, s_g, r_a, a, a_g, b, b_g, y_shock[indx],
                                        cost, util, deposit)

        cc  = kron(lambda, my_speye(I*J))
        ccu = kron(lambda, my_speye(I_g*J_g))

        @show cc[1:I*J+1,1:I*J+1]
@show lambda == Array(lambda')
@show lambda
        # full transition matrix
        A  = A + cc
        AT = (AT + ccu)'

        #----------------------------------------------------------------
        # KFE
        #----------------------------------------------------------------
        gIntermediate = dab_g_tilde_mat_inv * (AT * (dab_g_tilde_mat * gg)) + death_process * gg

        #----------------------------------------------------------------
        # Compute equilibrium conditions
        #----------------------------------------------------------------
        # HJB equation
        perm_mult   = (permanent==0) ? rrho + ddeath : rrho + ddeath - (1 - ggamma) * aggZ

        hjbResidual = vec(util.(c)) + A * vec(V) + V_Dot + VEErrors - perm_mult *
            reshape(V, I*J*N,1)

        # KFE
        gResidual = g_Dot - gIntermediate[1:n_g, 1]

        K_out = 0.0
        if K_liquid == 1
            K_out = sum((a_gg .+
                         b_gg) .* gg .* vec(dab_g))
        else
            K_out = sum(a_gg .* gg .* vec(dab_g))
        end

        K_Residual   = K_out - K
        r_b_out      = 0.0
        r_b_Residual = 0.0
        if r_b_fix      == 1
            r_b_out      = r_b_SS
            r_b_Residual = r_b_out - r_b
        elseif r_b_phi  == 1
            r_b_out      = sum(b_gg .* gg .* vec(dab_g))
            r_b_Residual = r_b_out - B_SS * exp(1/pphi * (r_b - r_b_SS))
        elseif B_fix    == 1
            # find death-corrected savings
            b_save       = dot(gIntermediate, dabb_g_aux)
            r_b_out      = 0.0
            r_b_Residual = r_b_out - b_save
        elseif K_liquid == 1
            r_b_out      = r_a_out - illiquid_wedge
            r_b_Residual = r_b_out - r_b
        end

        #if aggregate_variables == 1

            aggY_out = (K ^ aalpha) * (n_SS ^ (1 - aalpha))
            aggC_out = sum(vec(c_g) .* gg .* vec(dab_g))
        #=
        elseif distributional_variables == 1

            C_Var_out = sum(log(vec(c_g)).^2 .* gg .* vec(dab_g)) -
                sum(log(vec(c_g)) .* gg .* vec(dab_g)) ^ 2
            earn = log.((1-tau_I) * w * repeat(repeat(vec(y), inner=I_g), inner=J_g) .+
                         b_gg .* (repeat(repeat(r_b_g_vec, outer=J_g), outer=N) .+ ddeath*pam) .+
                         trans .+ a_gg .* (r_a + ddeath*pam))
            earn_Var_out = sum(vec(earn).^2 .* gg .* vec(dab_g)) -
                sum(vec(earn) .* gg .* vec(dab_g)) ^ 2

        elseif distributional_variables_1 == 1

            WHTM_indicator      = zeros(I_g,J_g,N)
            WHTM_indicator[b_g_0pos:b_g_0pos+1,a_g_0pos+2:end,:] .= 1.
            WHTM_out            = sum(vec(WHTM_indicator) .* gg .* vec(dab_g))
            C_WHTM_out          = sum(vec(WHTM_indicator) .* vec(c_g) .* gg .* vec(dab_g))

            PHTM_indicator      = zeros(I_g,J_g,N)
            PHTM_indicator[b_g_0pos:b_g_0pos+1,a_g_0pos:a_g_0pos+2:end,:] .= 1.
            PHTM_out            = sum(vec(PHTM_indicator) .* gg .* vec(dab_g))
            C_PHTM_out          = sum(vec(PHTM_indicator) .* vec(c_g) .* gg .* vec(dab_g))
        end
        =#
        Y_Residual        = 0.0
        C_Residual        = 0.0

        C_Var_Residual    = Array{Float64}(undef, 0)
        earn_Var_Residual = Array{Float64}(undef, 0)

        C_WHTM_Residual   = Array{Float64}(undef, 0)
        C_PHTM_Residual   = Array{Float64}(undef, 0)

        if aggregate_variables == 1
            Y_Residual        = aggY_out - aggY
            C_Residual        = aggC_out - aggC
        elseif distributional_variables == 1
            C_Var_Residual    = C_Var_out - C_Var
            earn_Var_Residual = earn_Var_out - earn_Var
        elseif distributional_variables_1 == 1
            C_WHTM_Residual   = C_WHTM_out - C_WHTM
            C_PHTM_Residual   = C_PHTM_out - C_PHTM
        end

        # Law of motion for aggregate tfp shock
        aggZ_Residual = aggZ_Dot - (-nnu_aggZ * aggZ + ssigma_aggZ * aggZ_Shock)

        # Return equilibrium conditions
        #if aggregate_variables == 1
            return [hjbResidual; gResidual; K_Residual; r_b_Residual; Y_Residual;
                    C_Residual; aggZ_Residual]
        #elseif distributional_variables == 1
        #    return [hjbResidual; gResidual; K_Residual; r_b_Residual; C_Var_Residual;
        #            earn_Var_Residual; aggZ_Residual]
        #elseif distributional_variables_1 == 1
        #    return [hjbResidual; gResidual; K_Residual; r_b_Residual; C_WHTM_Residual;
        #            C_PHTM_Residual; aggZ_Residual]
        #end
        #return [hjbResidual; gResidual; K_Residual; r_b_Residual; aggZ_Residual]
    end

    out = get_residuals_lite(zeros(Float64, 2 * nVars + nEErrors + 1))
    #JLD2.jldopen("/home/rcerxs30/.julia/dev/DSGE/src/models/heterogeneous_agent/two_asset_hank/eqcond_after_.jld2", true, true, true, IOStream) do file
    #    file["residuals"] = out
    #end
    test_out = load("/home/rcerxs30/.julia/dev/DSGE/src/models/heterogeneous_agent/two_asset_hank/eqcond_after_1e12.jld2", "residuals")
    @assert test_out == out
    @time get_residuals_lite(zeros(Float64, 2 * nVars + nEErrors + 1))

    x = zeros(Float64, 2 * nVars + nEErrors + 1)
    @time derivs = ForwardDiff.sparse_jacobian(get_residuals_lite, x)

    nstates = nVars # n_states(m)
    n_s_exp = nEErrors # n_shocks_expectational(m)
    n_s_exo = n_Z # n_shocks_exogenous(m)
    #vars = zeros(Float64, 2 * nstates + n_s_exp + n_s_exo)

    Γ1 = -derivs[:, 1:nstates]
    Γ0 =  derivs[:,   nstates +           1:2 * nstates]
    Π  = -derivs[:, 2*nstates +           1:2 * nstates + n_s_exp]
    Ψ  = -derivs[:, 2*nstates + n_s_exp + 1:2 * nstates + n_s_exp + n_s_exo]
    C  = spzeros(Float64, nstates)

    if typeof(Ψ) == Vector{Float64}
        Ψ = reshape(Ψ, length(Ψ), n_s_exo)
    end

    test_out = load("/home/rcerxs30/.julia/dev/DSGE/src/models/heterogeneous_agent/two_asset_hank/data/eqcond_output_matlab.jld2")
    @show test_out["g0"] == Γ0
    @show isapprox(test_out["g1"], Γ1, rtol=1e-4)
    @show test_out["psi"] == Ψ
    @show test_out["pi"] == Π
    @show test_out["constant"] == C

    return Γ0, Γ1, Ψ, Π, C
end
