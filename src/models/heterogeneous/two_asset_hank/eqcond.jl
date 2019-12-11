using Statistics, SparseArrays
using MAT, DelimitedFiles
"""
```
eqcond(m::TwoAssetHANK)
```
Expresses the equilibrium conditions in canonical form using Γ0, Γ1, C, Ψ, and Π matrices.

# Outputs
* `Γ0` (`n_states` x `n_states`) holds coefficients of current time states.
* `Γ1` (`n_states` x `n_states`) holds coefficients of lagged states.
* `C`  (`n_states` x `1`) is a vector of constants
* `Ψ`  (`n_states` x `n_shocks_exogenous`) holds coefficients of iid shocks.
* `Π`  (`n_states` x `n_states_expectational`) holds coefficients of expectational states.
"""
function eqcond(m::TwoAssetHANK)

    # Read in parameters
    α  = m[:aalpha].value::Float64
    δ  = m[:ddelta].value::Float64
    ρ  = m[:rrho].value::Float64
    γ  = m[:ggamma].value::Float64
    κ  = m[:kappa].value::Float64
    ξ  = m[:xxi].value::Float64
    ϕ  = m[:pphi].value::Float64
    χ0 = m[:chi0].value::Float64
    χ1 = m[:chi1].value::Float64
    χ2 = m[:chi2].value::Float64

    ddeath = m[:ddeath].value::Float64
    a_lb   = m[:a_lb].value::Float64
    pam    = m[:pam].value::Float64
    τ_I    = m[:tau_I].value::Float64
    trans  = m[:trans].value::Float64
    n_SS   = m[:n_SS].value::Float64
    ν_aggZ = m[:nnu_aggZ].value::Float64
    σ_aggZ = m[:ssigma_aggZ].value::Float64

    # Set liquid rates
    r_b_SS       = m[:r_b_SS].value::Float64
    r_b_borr_SS  = m[:r_b_borr_SS].value::Float64
    borrwedge_SS = m[:borrwedge_SS].value::Float64

    K_liquid = get_setting(m, :K_liquid)::Bool
    r_b_fix  = get_setting(m, :r_b_fix)::Bool
    r_b_phi  = get_setting(m, :r_b_phi)::Bool
    B_fix    = get_setting(m, :B_fix)::Bool

    aggregate_variables        = get_setting(m, :aggregate_variables)::Bool
    distributional_variables   = get_setting(m, :distributional_variables)::Bool
    distributional_variables_1 = get_setting(m, :distributional_variables_1)::Bool
    permanent                  = get_setting(m, :permanent)::Bool

    I   = get_setting(m, :I)::Int64
    J   = get_setting(m, :J)::Int64
    I_g = get_setting(m, :I_g)::Int64
    J_g = get_setting(m, :J_g)::Int64
    N   = get_setting(m, :y_size)::Int64

    a   = get_setting(m, :a)::Vector{Float64}
    b   = get_setting(m, :b)::Vector{Float64}
    a_g = get_setting(m, :a_g)::Vector{Float64}
    b_g = get_setting(m, :b_g)::Vector{Float64}

    agrid_new = get_setting(m, :agrid_new)::Int64
    bgrid_new = get_setting(m, :bgrid_new)::Int64
    ygrid_new = get_setting(m, :ygrid_new)::Int64

    amin = get_setting(m, :amin)
    bmin = get_setting(m, :bmin)
    amax = get_setting(m, :amax)
    bmax = get_setting(m, :bmax)

    KL      = get_setting(m, :KL_0)::Float64

    n_v = get_setting(m, :n_v)::Int64
    n_g = get_setting(m, :n_g)::Int64
    n_p = get_setting(m, :n_p)::Int64
    n_Z = get_setting(m, :n_Z)::Int64
    nVars    = get_setting(m, :nVars)::Int64
    nEErrors = get_setting(m, :nEErrors)::Int64

    #--- Taken from what used to be inside
    a, a_g, a_g_0pos_arr          = create_a_grid(agrid_new, J, J_g, amin, amax)
    _, _, _, b, b_g, b_g_0pos_arr = create_b_grid(bgrid_new, I, I_g)
    lambda, y, y_mean, y_dist, _  = create_y_grid(N, ygrid_new)
    a_g_0pos = first(a_g_0pos_arr)
    b_g_0pos = first(b_g_0pos_arr)

    # Construct problem functions
    util, deposit, cost = construct_problem_functions(γ, χ0, χ1, χ2, a_lb)

    vars_SS     = vec(m[:vars_SS].value)::Vector{Float64}
    #vars_SS = vec(load("data/vars_SS.jld2","vars_SS"))
    V_SS   = vars_SS[1:n_v]
    g_SS   = vars_SS[n_v + 1 : n_v + n_g]
    K_SS   = vars_SS[n_v + n_g + 1]
    r_b_SS = vars_SS[n_v + n_g + 2]

    # Aggregate output and aggregate consumption
    aggY_SS = aggregate_variables ? vars_SS[n_v+n_g+3] : 0.0
    aggC_SS = aggregate_variables ? vars_SS[n_v+n_g+4] : 0.0

    # Consumption and earnings inequality
    C_Var_SS    = distributional_variables ? vars_SS[n_v+n_g+3] : 0.0
    earn_Var_SS = distributional_variables ? vars_SS[n_v+n_g+4] : 0.0

    # Consumption of wealthy and poor hand-to-mouth
    C_WHTM_SS = distributional_variables_1 ? vars_SS[n_v+n_g+3] : 0.0
    C_PHTM_SS = distributional_variables_1 ? vars_SS[n_v+n_g+4] : 0.0

    aggZ_SS = vars_SS[n_v+n_g+n_p+1] # Aggregate Z

    @inline function get_residuals(vars::Vector{T}) where {T<:Real}

        # ------- Unpack variables -------
        V      = reshape(vars[1:n_v] .+ V_SS, I, J, N) # value function
        g      = vars[n_v + 1 : n_v + n_g] .+ g_SS     # distribution
        K::T   = vars[n_v + n_g + 1]        + K_SS     # aggregate capital
        r_b::T = vars[n_v + n_g + 2]        + r_b_SS

        # Aggregate output, aggregate consumption
        aggY::T = aggregate_variables ? vars[n_v+n_g+3] + aggY_SS : zero(T)
        aggC::T = aggregate_variables ? vars[n_v+n_g+4] + aggC_SS : zero(T)

        # Consumption and earnings inequality
        C_Var::T    = distributional_variables ? vars[n_v+n_g+3] + C_Var_SS    : zero(T)
        earn_Var::T = distributional_variables ? vars[n_v+n_g+4] + earn_Var_SS : zero(T)

        # Consumption of wealthy and poor hand-to-mouth
        C_WHTM::T = distributional_variables_1 ? vars[n_v+n_g+3] + C_WHTM_SS : zero(T)
        C_PHTM::T = distributional_variables_1 ? vars[n_v+n_g+4] + C_PHTM_SS : zero(T)

        aggZ::T       = vars[n_v + n_g + n_p + 1] + aggZ_SS
        V_Dot         = vars[nVars +       1 : nVars + n_v]
        g_Dot         = vars[nVars + n_v + 1 : nVars + n_v + n_g]
        aggZ_Dot::T   = vars[nVars + n_v + n_g + n_p + 1]
        VEErrors      = vars[2*nVars + 1 : 2 * nVars + n_v]
        aggZ_Shock::T = vars[2*nVars + nEErrors + 1]

        # Prices
        w   = (1-α) * (K ^ α) * n_SS ^ (-α) * (!permanent ? exp(aggZ) ^ (1-α) : 1.0)
        r_a = α * (K ^ (α-1)) * ((!permanent ? exp(aggZ) : 1.0) * n_SS) ^ (1-α) - δ

        # Set liquid rates
        r_b_borr = r_b .+ borrwedge_SS

        # Initialize grids
        dab_g = reshape(repeat(backward_difference(a_g, b_g), N, 1), I_g, J_g, N)
        a_gg  = repeat(repeat(a_g, inner=I_g), outer=N)
        b_gg  = repeat(repeat(b_g, outer=J_g), outer=N)

        g_end = (1 - sum(g .* vec(dab_g)[1:end-1])) / dab_g[I_g, J_g, N]
        gg    = vcat(g, g_end)

        loc = findall(b .== 0)
        dab_g_small          = vec(dab_g[:,:,1])
        dab_g_small          = dab_g_small ./ dab_g_small[loc] * ddeath
        dab_g_small[loc]    .= 0.0
        death_process        = -ddeath * my_speye(I_g * J_g)
        death_process[loc,:] = vec(dab_g_small)
        death_process        = kron(my_speye(N), death_process)

        r_b_vec = r_b .* (b .>= 0) + r_b_borr .* (b .< 0)

        # Other necessary objects
        y_shock      = y .* exp.(κ * aggZ * (y .- y_mean) ./ std(y))
        y_shock_mean = dot(y_shock, y_dist)
        y_shock      = real(y_shock ./ y_shock_mean .* y_mean)

        println("Timing: solve_hjb()")
        @time c, s, d  = solve_hjb(V, a_lb, γ, ddeath, pam, trans, ξ, τ_I, aggZ, w,
                                   r_b_vec, y_shock, a, b, cost, util, deposit;
                                   permanent=permanent)

        interp_decision = kron(my_speye(N), interp(b_g, a_g, b, a))
        d_g = reshape(interp_decision * vec(d), I_g, J_g, N)
        s_g = reshape(interp_decision * vec(s), I_g, J_g, N)
        c_g = reshape(interp_decision * vec(c), I_g, J_g, N)

        # Derive transition matrices
        println("Timing: transition()")
        @time A, AT = transition(ddeath, pam, ξ, w, a_lb, aggZ, d, d_g, s, s_g, r_a,
                                 a, a_g, b, b_g, y_shock, cost; permanent=permanent)

        cc  = kron(lambda, my_speye(I*J))
        ccu = kron(lambda, my_speye(I_g*J_g))

        # Full transition matrix
        A  = A + cc
        AT = (AT + ccu)'

        #----------------------------------------------------------------
        # KFE
        #----------------------------------------------------------------
        ggdab_g = gg .* vec(dab_g)
        gIntermediate = (1.0 ./ vec(dab_g)) .* (AT * ggdab_g) + death_process * gg

        #----------------------------------------------------------------
        # Compute equilibrium conditions
        #----------------------------------------------------------------
        # HJB equation
        perm_mult   = !permanent ? ρ + ddeath : ρ + ddeath - (1 - γ) * aggZ
        hjbResidual = vec(util.(c)) + A * vec(V) + V_Dot + VEErrors - perm_mult * vec(V)

        # KFE
        gResidual  = g_Dot - gIntermediate[1:n_g]
        K_Residual = (K_liquid ? sum((a_gg .+ b_gg) .* ggdab_g) : sum(a_gg .* ggdab_g)) - K
        r_b_Residual = if r_b_fix
                           r_b_SS - r_b
                       elseif r_b_phi
                           sum(b_gg .* ggdab_g) - B_SS * exp(1/ϕ * (r_b - r_b_SS))
                       elseif B_fix
                           -dot(gIntermediate, vec(dab_g) .* b_gg)
                       elseif K_liquid
                           r_a_out - illiquid_wedge - r_b
                       end

        # Law of motion for aggregate TFP shock
        aggZ_Residual = aggZ_Dot - (-ν_aggZ * aggZ + σ_aggZ * aggZ_Shock)

        # Return equilibrium conditions
        if aggregate_variables

            aggY_out = (K ^ α) * (!permanent ? exp(aggZ) * n_SS : n_SS ) ^ (1 - α)
            aggC_out = sum(vec(c_g) .* ggdab_g)
            Y_Residual = aggY_out - aggY
            C_Residual = aggC_out - aggC

            return [hjbResidual; gResidual; K_Residual; r_b_Residual; Y_Residual;
                    C_Residual; aggZ_Residual]

        elseif distributional_variables

            r_b_g = repeat(repeat(r_b .* (b_g .>= 0) + r_b_borr .* (b_g .< 0), outer=J_g), outer=N)
            C_Var_out = sum(log.(vec(c_g)).^2 .* ggdab_g) - sum(log.(vec(c_g)) .* ggdab_g) ^ 2

            earn = vec(log.((1-τ_I) * w * repeat(repeat(vec(y), inner=I_g), inner=J_g) .+ b_gg .*
                            r_b_g .+ ddeath*pam) .+ trans .+ a_gg .* (r_a + ddeath*pam))
            earn_Var_out = sum(earn.^2 .* ggdab_g) - sum(earn .* ggdab_g) ^ 2

            C_Var_Residual    = C_Var_out - C_Var
            earn_Var_Residual = earn_Var_out - earn_Var

            return [hjbResidual; gResidual; K_Residual; r_b_Residual; C_Var_Residual;
                    earn_Var_Residual; aggZ_Residual]

        elseif distributional_variables_1

            WHTM_indicator  = zeros(I_g,J_g,N)
            WHTM_indicator[b_g_0pos:b_g_0pos+1, a_g_0pos+2:end, :] .= 1.
            C_WHTM_out      = sum(vec(WHTM_indicator) .* vec(c_g) .* ggdab_g)

            PHTM_indicator  = zeros(I_g,J_g,N)
            PHTM_indicator[b_g_0pos:b_g_0pos+1, a_g_0pos:a_g_0pos+2:end, :] .= 1.
            C_PHTM_out      = sum(vec(PHTM_indicator) .* vec(c_g) .* ggdab_g)

            C_WHTM_Residual = C_WHTM_out - C_WHTM
            C_PHTM_Residual = C_PHTM_out - C_PHTM

            return [hjbResidual; gResidual; K_Residual; r_b_Residual; C_WHTM_Residual;
                    C_PHTM_Residual; aggZ_Residual]
        end
        return [hjbResidual; gResidual; K_Residual; r_b_Residual; aggZ_Residual]
    end

    out = get_residuals(zeros(Float64, 2 * nVars + nEErrors + 1))

    #JLD2.jldopen("reference/eqcond_after_.jld2", true, true, true, IOStream) do file
    #    file["residuals"] = out
    #end

    my_out = vec(DelimitedFiles.readdlm("reference/my_residuals.csv", ','))

    @show maximum(abs.(my_out - vec(out)))#, length(findall(x->x>1e-5, abs.(my_out - out)))
    #@assert isapprox(my_out, vec(out), rtol=1e-4)

    nstates = nVars    # n_states(m)
    n_s_exp = nEErrors # n_shocks_expectational(m)
    n_s_exo = n_Z      # n_shocks_exogenous(m)

    @time get_residuals(zeros(Float64, 2 * nstates + n_s_exp + 1))

    x = zeros(Float64, 2 * nstates + n_s_exp + 1)
    @time derivs = ForwardDiff.sparse_jacobian(get_residuals, x)

    # vars = zeros(Float64, 2 * nstates + n_s_exp + n_s_exo)

    Γ1 = -derivs[:, 1:nstates]
    Γ0 =  derivs[:,   nstates +           1:2 * nstates]
    Π  = -derivs[:, 2*nstates +           1:2 * nstates + n_s_exp]
    Ψ  = -derivs[:, 2*nstates + n_s_exp + 1:2 * nstates + n_s_exp + n_s_exo]
    C  =  spzeros(Float64, nstates)

    if typeof(Ψ) == Vector{Float64}
        Ψ = reshape(Ψ, length(Ψ), n_s_exo)
    end

    test_out = load("data/eqcond_output_matlab.jld2")
    @assert test_out["g0"] == Γ0
    @show maximum(abs.(test_out["g1"] - Γ1))
    @assert test_out["g1"] ≈ Γ1
    @assert test_out["psi"] == Ψ
    @assert test_out["pi"] == Π

    return Γ0, Γ1, Ψ, Π, C
end
