using Statistics, SparseArrays
using MAT, DelimitedFiles, ForwardDiff
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
    r_b_borr_SS  = m[:r_b_borr_SS].value::Float64
    borrwedge_SS = m[:borrwedge_SS].value::Float64

    lambda::Matrix{Float64} = get_setting(m, :lambda)

    K_liquid                   = get_setting(m, :K_liquid)::Bool
    aggregate_variables        = get_setting(m, :aggregate_variables)==1#::Bool
    distributional_variables   = get_setting(m, :distributional_variables)==1#::Int64
    distributional_variables_1 = get_setting(m, :distributional_variables_1)==1#::Int64
    permanent                  = get_setting(m, :permanent)::Bool

    I      = get_setting(m, :I)::Int64
    J      = get_setting(m, :J)::Int64
    I_g    = get_setting(m, :I_g)::Int64
    J_g    = get_setting(m, :J_g)::Int64
    N      = get_setting(m, :N)::Int64

    a      = get_setting(m, :a)::Vector{Float64}
    b      = get_setting(m, :b)::Vector{Float64}
    a_g    = get_setting(m, :a_g)::Vector{Float64}
    b_g    = get_setting(m, :b_g)::Vector{Float64}

    agrid_new = get_setting(m, :agrid_new)::Int64
    bgrid_new = get_setting(m, :bgrid_new)::Int64
    ygrid_new = get_setting(m, :ygrid_new)::Int64

    amin = get_setting(m, :amin)
    bmin = get_setting(m, :bmin)
    amax = get_setting(m, :amax)
    bmax = get_setting(m, :bmax)

    y      = vec(get_setting(m, :y))::Vector{Complex{Float64}}
    y_dist = get_setting(m, :y_dist)::Vector{Complex{Float64}}
    y_mean = get_setting(m, :y_mean)::Complex{Float64}
    KL      = get_setting(m, :KL_0)::Float64
    r_b_fix = get_setting(m, :r_b_fix) ? 1 : 0::Int64

    #--- Taken from what used to be inside
    #a, a_g, a_g_0pos          = create_a_grid(agrid_new, J, J_g, amin, amax)
    #    _, _, _, b, b_g, b_g_0pos = create_b_grid(bgrid_new, I, I_g)
    #lambda, y, y_mean, y_dist, _ = create_y_grid(N, ygrid_new)

    n_v = get_setting(m, :n_v)::Int64
    n_g = get_setting(m, :n_g)::Int64
    n_p = get_setting(m, :n_p)::Int64
    n_Z = get_setting(m, :n_Z)::Int64
    nVars    = get_setting(m, :nVars)::Int64
    nEErrors = get_setting(m, :nEErrors)::Int64

    #vars_SS     = vec(m[:vars_SS].value)::Vector{Float64}
    vars_SS = vec(load("data/vars_SS.jld2","vars_SS"))
    V_SS   = reshape(vars_SS[1:n_v], I*J, N)
    g_SS   = vars_SS[n_v + 1 : n_v + n_g]
    K_SS   = vars_SS[n_v + n_g + 1]
    r_b_SS = vars_SS[n_v + n_g + 2]

    # Aggregate output and aggregate consumption
    aggY_SS = aggregate_variables ? vars_SS[n_v+n_g+3] : 0.0 # aggregate output
    aggC_SS = aggregate_variables ? vars_SS[n_v+n_g+4] : 0.0 # aggregate consumption

    # Consumption and earnings inequality
    C_Var_SS    = distributional_variables ? vars_SS[n_v+n_g+3] : 0.0
    earn_Var_SS = distributional_variables ? vars_SS[n_v+n_g+4] : 0.0

    # Consumption of wealthy and poor hand-to-mouth
    C_WHTM_SS = distributional_variables_1 ? vars_SS[n_v+n_g+3] : 0.0
    C_PHTM_SS = distributional_variables_1 ? vars_SS[n_v+n_g+4] : 0.0

    aggZ_SS = vars_SS[n_v+n_g+n_p+1] # Aggregate Z

    # Construct problem functions
    util, deposit, cost = construct_problem_functions(ggamma, chi0, chi1, chi2, a_lb)

    @inline function get_residuals_lite(vars::Vector{T}) where {T<:Real}
        # ------- Unpack variables -------
        V::Matrix{T} = reshape(vars[1:n_v] .+ V_SS[:,indx], I, J)  # Value function
        g::Vector{T} = vars[n_v + 1 : n_v + n_g] .+ g_SS           # Distribution
        K::T         = vars[n_v + n_g + 1] + K_SS                  # Aggregate capital
        r_b::T       = vars[n_v + n_g + 2] + r_b_SS

        # Aggregate output, aggregate consumption
        aggY::T     = aggregate_variables      ? vars[n_v+n_g+3] + aggY_SS     : zero(T)
        aggC::T     = aggregate_variables      ? vars[n_v+n_g+4] + aggC_SS     : zero(T)
        C_Var::T    = distributional_variables ? vars[n_v+n_g+3] + C_Var_SS    : zero(T)
        earn_Var::T = distributional_variables ? vars[n_v+n_g+4] + earn_Var_SS : zero(T)
        C_WHTM::T = distributional_variables_1 ? vars[n_v+n_g+3] + C_WHTM_SS   : zero(T)
        C_PHTM::T = distributional_variables_1 ? vars[n_v+n_g+4] + C_PHTM_SS   : zero(T)

        # Aggregate Z
        aggZ::T             = vars[n_v+n_g+n_p+1] + aggZ_SS
        V_Dot::Vector{T}    = vars[    nVarsi +       1 : nVarsi + n_v]
        g_Dot::Vector{T}    = vars[    nVarsi + n_v + 1 : nVarsi + n_v + n_g]
        aggZ_Dot::T         = vars[    nVarsi + n_v + n_g + n_p + 1]
        VEErrors::Vector{T} = vars[2 * nVarsi + 1 : 2 * nVarsi + n_v]
        aggZ_Shock::T       = vars[2 * nVarsi + nEErrors + 1]

        # Prices, liquid rates
        w   = (1 - aalpha) * (K ^ aalpha) * n_SS ^ (-aalpha) *
            (!permanent ? exp(aggZ) ^ (1-aalpha) : 1.)
        r_a = aalpha * (K ^ (aalpha - 1)) * ((!permanent ? exp(aggZ) : 1.0) *
            n_SS) ^ (1 - aalpha) - ddelta
        r_b_borr = r_b .+ borrwedge_SS

        # Set grid
        dab_g_tilde = backward_difference(a_g, b_g)

        dab_g_o = reshape(repeat(reshape(dab_g_tilde, I_g, J_g), N, 1), I_g, J_g, N)
        g_end = (1 - sum(g .* vec(dab_g_o)[1:end-1])) / dab_g_tilde[I_g*J_g]
        @show size(dab_g_o)
        @show size(g .* vec(dab_g_o)[1:end-1])
        @show size(vec(broadcast(*, dab_g_tilde, reshape(vcat(g, 0.0), I_g*J_g, N)))[1:end-1])

        g      = vcat(g, g_end)
        g_inds = reshape(g, I_g*J_g, N)
        gg     = g_inds[:,indx]

        # Other necessary objects
        y_shock      = y .* exp.(kappa * aggZ * (y .- y_mean) ./ std(y))
        y_shock_mean = dot(y_shock, y_dist)
        y_shock      = real(y_shock ./ y_shock_mean .* y_mean)

        println("Timing: solve_hjb()")
        @time c, s, d  = solve_hjb_lite(V, a_lb, ggamma, permanent, ddeath, pam, aggZ,
                                        xxi, tau_I, w, trans, r_b, r_b_borr, y_shock[indx],
                                        a, b, cost, util, deposit)

        @time interp_decision = interp(b_g, a_g, b, a)
        d_g = reshape(interp_decision * vec(d), I_g, J_g)
        s_g = reshape(interp_decision * vec(s), I_g, J_g)
        c_g = reshape(interp_decision * vec(c), I_g, J_g)

        @show size(interp_decision)
        # Derive transition matrices
        println("Timing: transition_deriva()")
        @time A, AT = transition_deriva_lite(permanent, ddeath, pam, xxi, w, a_lb, aggZ,
                                        d, d_g, s, s_g, r_a, a, a_g, b, b_g, y_shock[indx],
                                        lambda[indx,indx], cost, util, deposit)

        #----------------------------------------------------------------
        # KFE
        #----------------------------------------------------------------
        @time gIntermediate, cc = solve_kfe(dab_g_tilde, b, g_inds, gg, ddeath,
                                            lambda, AT, indx, V_SS)

        #----------------------------------------------------------------
        # Compute equilibrium conditions
        #----------------------------------------------------------------
        # HJB equation
        perm_mult   = !permanent ? rrho + ddeath : rrho + ddeath - (1 - ggamma) * aggZ
        hjbResidual = vec(util.(c)) + A * vec(V) + cc + V_Dot + VEErrors - perm_mult * vec(V)

        # KFE
        a_gg = repeat(a_g, inner=I_g)
        b_gg = repeat(b_g, outer=J_g)

        @show length(g_Dot), length(gIntermediate)
        gResidual  = g_Dot[I_g*J_g*(indx-1)+1 : I_g*J_g*indx] - gIntermediate
        K_Residual = (K_liquid ? sum((a_gg .+ b_gg) .* gg .* dab_g_tilde) :
                                 sum( a_gg          .* gg .* dab_g_tilde)) - K
        r_b_Residual = 0.0
        if r_b_fix      == 1
            r_b_out      = r_b_SS
            r_b_Residual = r_b_out - r_b
        elseif r_b_phi  == 1
            r_b_out      = sum(b_gg .* gg .* dab_g_tilde)
            r_b_Residual = r_b_out - B_SS * exp(1 / pphi * (r_b - r_b_SS))
        elseif B_fix    == 1
            # Find death-corrected savings
            b_save       = dot(gIntermediate, dab_g_tilde .* b_gg)
            r_b_out      = 0.0
            r_b_Residual = r_b_out - b_save
        elseif K_liquid
            r_b_out      = r_a_out - illiquid_wedge
            r_b_Residual = r_b_out - r_b
        end

        Y_Residual        = zero(T)
        C_Residual        = zero(T)

        C_Var_Residual    = Array{T}(undef, 0)
        earn_Var_Residual = Array{T}(undef, 0)

        C_WHTM_Residual   = Array{T}(undef, 0)
        C_PHTM_Residual   = Array{T}(undef, 0)

        if aggregate_variables

            aggY_out = (K ^ aalpha) * (!permanent ? exp(aggZ) * n_SS : n_SS ) ^ (1 - aalpha)
            aggC_out = sum(vec(c_g) .* gg .* dab_g_tilde)
            Y_Residual = aggY_out - aggY
            C_Residual = aggC_out - aggC

        elseif distributional_variables
            r_b_g_vec = r_b .* (b_g .>= 0) + r_b_borr .* (b_g .< 0)
            C_Var_out = sum(log(vec(c_g)).^2 .* gg .* dab_g_tilde) -
                sum(log(vec(c_g)) .* gg .* dab_g_tilde) ^ 2
            earn = log.((1-tau_I) * w * repeat(repeat(vec(y), inner=I_g), inner=J_g) .+
                         b_gg .* (repeat(repeat(r_b_g_vec, outer=J_g), outer=N) .+ ddeath*pam) .+
                         trans .+ a_gg .* (r_a + ddeath*pam))
            earn_Var_out = sum(vec(earn).^2 .* gg .* dab_g_tilde) -
                sum(vec(earn) .* gg .* dab_g_tilde) ^ 2

            C_Var_Residual    = C_Var_out - C_Var
            earn_Var_Residual = earn_Var_out - earn_Var

        elseif distributional_variables_1

            WHTM_indicator      = zeros(I_g,J_g,N)
            WHTM_indicator[b_g_0pos:b_g_0pos+1,a_g_0pos+2:end,:] .= 1.
            WHTM_out            = sum(vec(WHTM_indicator) .* gg .* dab_g_tilde)
            C_WHTM_out          = sum(vec(WHTM_indicator) .* vec(c_g) .* gg .* dab_g_tilde)

            PHTM_indicator      = zeros(I_g,J_g,N)
            PHTM_indicator[b_g_0pos:b_g_0pos+1,a_g_0pos:a_g_0pos+2:end,:] .= 1.
            PHTM_out            = sum(vec(PHTM_indicator) .* gg .* dab_g_tilde)
            C_PHTM_out          = sum(vec(PHTM_indicator) .* vec(c_g) .* gg .* dab_g_tilde)

            C_WHTM_Residual = C_WHTM_out - C_WHTM
            C_PHTM_Residual = C_PHTM_out - C_PHTM
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

    indx = 1
    n_v = Int(n_v/N)
    #n_g = Int((n_g + 1)/N)
    nVarsi = n_v + n_g + n_p + 1

    out = get_residuals_lite(zeros(Float64, 2 * nVarsi + nEErrors + 1))
    @time out = get_residuals_lite(zeros(Float64, 2 * nVarsi + nEErrors + 1))

    #JLD2.jldopen("reference/eqcond_after_.jld2", true, true, true, IOStream) do file
    #    file["residuals"] = out
    #end

    my_out = vec(DelimitedFiles.readdlm("reference/my_residuals.csv", ','))

    @assert maximum(abs.(my_out[1:2000] - vec(out)[1:2000])) < 1e-7
    @show maximum(abs.(my_out[1:2000] - vec(out)[1:2000]))
    @assert isapprox(my_out[60001:62250], vec(out)[2001:4250], rtol=1e-4)

    x = zeros(Float64, 2 * nVarsi + nEErrors + 1)
    @time derivs = ForwardDiff.sparse_jacobian(get_residuals_lite, x)

    nstates = nVars    # n_states(m)
    n_s_exp = nEErrors # n_shocks_expectational(m)
    n_s_exo = n_Z      # n_shocks_exogenous(m)

    # vars = zeros(Float64, 2 * nstates + n_s_exp + n_s_exo)

    Γ1 = spzeros(Float64, nstates, nstates)#-derivs[:, 1:nstates]
    Γ0 = spzeros(Float64, nstates, nstates)# derivs[:,   nstates +  1:2 * nstates]
    Π  = spzeros(Float64, nstates, n_s_exp)#-derivs[:, 2*nstates +  1:2 * nstates + n_s_exp]
    Ψ  = spzeros(Float64, nstates, n_s_exo)#-derivs[:, 2*nstates + n_s_exp + 1:2 * nstates + n_s_exp + n_s_exo]
    C  =  spzeros(Float64, nstates)

    for i=1:1
        indx = i
        x = zeros(Float64, 2 * nVarsi + nEErrors + 1)
        @time derivs = ForwardDiff.sparse_jacobian(get_residuals_lite, x)
        #Γ1[] = -
        #Γ0[] =
        #Π[]  =
        #Ψ[]  =
    end

    if typeof(Ψ) == Vector{Float64}
        Ψ = reshape(Ψ, length(Ψ), n_s_exo)
    end

    test_out = load("data/eqcond_output_matlab.jld2")
    @assert test_out["g0"] == Γ0
    @assert test_out["g1"] ≈ Γ1
    @assert test_out["psi"] == Ψ
    @assert test_out["pi"] == Π

    return Γ0, Γ1, Ψ, Π, C
end
