# Set up model in canonical form
function equilibrium_conditions(vars_SS::OrderedDict{Symbol, Any},
                                grids::Dict{Symbol, Any},
                                params::Dict{Symbol, Float64},
                                approx_params::Dict{Symbol, Any})

    n_vars = grids[:n_vars]
    n_exp_errors = grids[:n_exp_errors]
    n_shocks = grids[:n_shocks]

    x = zeros(Float64, 2*n_vars + n_exp_errors + n_shocks)

    function f{T<:Real}(x::Vector{T}; vars_ss = deepcopy(vars_SS))
        v_residual = get_residuals(x, vars_ss, grids, params, approx_params)
        return vcat(values(v_residual)...)
    end

    derivs = ForwardDiff.jacobian(f, x)
    Γ1 = -derivs[:,1:n_vars];
    Γ0 = derivs[:,n_vars+1:2*n_vars];
    Π = -derivs[:,2*n_vars+1:2*n_vars+n_exp_errors];
    Ψ = -derivs[:,2*n_vars+n_exp_errors+1:2*n_vars+n_exp_errors+n_shocks];
    C = zeros(n_vars, 1)

    canonical_form = Dict{Symbol, Any}()
    canonical_form[:Γ1] = Γ1
    canonical_form[:Γ0] = Γ0
    canonical_form[:Π] = Π
    canonical_form[:Ψ] = Ψ
    canonical_form[:C] = C

    return canonical_form
end

function get_residuals(vars_SS::OrderedDict{Symbol, Any},
                       grids::Dict{Symbol, Any},
                       params::Dict{Symbol, Float64},
                       approx_params::Dict{Symbol, Any})

    n_vars = grids[:n_vars]
    n_exp_errors = grids[:n_exp_errors]
    n_shocks = grids[:n_shocks]
    x = zeros(Float64, 2*n_vars + n_exp_errors + n_shocks)

    return get_residuals(x, vars_SS, grids, params, approx_params)
end

function get_residuals{T<:Real}(x::Vector{T},
                                vars_SS::OrderedDict{Symbol, Any},
                                grids::Dict{Symbol, Any},
                                params::Dict{Symbol, Float64},
                                approx_params::Dict{Symbol, Any})

    # Read in steady state vars, params, and grids
    niter_hours = approx_params[:niter_hours]

    n_v = grids[:n_jump_vars]
    n_g = grids[:n_state_vars]
    n_p = grids[:n_static_conditions]
    n_shocks = grids[:n_shocks]
    n_exp_errors = grids[:n_exp_errors]
    n_vars = grids[:n_vars]

    I, J = size(grids[:zz])

    vars_SS[:V_SS] = vars_SS[:V_SS] + reshape(x[1:n_v - 1], I, J)
    vars_SS[:inflation_SS] = vars_SS[:inflation_SS] + x[n_v]

    gg = vec(vars_SS[:g_SS])[1:end-1] + x[n_v + 1 : n_v + n_g - 1]
    vars_SS[:g_SS] = reshape(vcat(gg, vars_SS[:g_SS][end,end]), I, J)

    vars_SS[:w_SS] = vars_SS[:w_SS] + x[n_v + n_g + 1]
    vars_SS[:N_SS] = vars_SS[:N_SS] + x[n_v + n_g + 2]
    vars_SS[:C_SS] = vars_SS[:C_SS] + x[n_v + n_g + 3]
    vars_SS[:Y_SS] = vars_SS[:Y_SS] + x[n_v + n_g + 4]
    vars_SS[:B_SS] = vars_SS[:B_SS] + x[n_v + n_g + 5]

    V  = vars_SS[:V_SS]
    w = vars_SS[:w_SS]
    h_SS = vars_SS[:h_SS]
    gg = vec(vars_SS[:g_SS])[1:end-1]

    a   = grids[:a]
    zz  = grids[:zz]
    ymarkov_combined = grids[:ymarkov_combined]
    amax = maximum(a)
    amin = minimum(a)

    coefrra = params[:coefrra]
    maxhours   = params[:maxhours]

    dazf, dazb, azdelta, aa, adelta, azdelta, azdelta_mat = initialize_diff_grids(zz, a)
    azdeltavec = vec(azdelta)

    g_end = (1 - sum(gg .* azdeltavec[1:end-1])) ./ azdeltavec[end]
    g = vcat(gg, g_end)
    MP = 0. + x[n_v + n_g] # @

    VDot = zeros(n_v - 1) + x[n_vars + 1 : n_vars + n_v - 1] # @
    inflationDot = 0. + x[n_vars + n_v] # @
    gDot = zeros(n_g - 1)  + x[n_vars + n_v + 1 : n_vars + n_v + n_g - 1] # @
    MPDot = 0. + x[n_vars + n_g + n_v] # @

    VEErrors = zeros(n_v - 1) + x[2*n_vars + 1 : 2*n_vars + n_v - 1] # @
    inflationEError = 0. + x[2*n_vars + n_v] # @

    MPShock = 0. + x[2*n_vars + n_v + 1] # @

    TFP = 1.

    # Initialize other variables, using V to ensure everything is a dual number
    # Note: Ensure that copying a dual number still gives a dual number in memory
    Vaf = copy(V)
    Vab = copy(V)
    cf  = similar(Vaf)
    hf  = similar(Vaf)
    sf  = similar(Vaf)

    cb  = similar(Vab)
    hb  = similar(Vab)
    sb  = similar(Vab)

    c0 = similar(cb)

    A = zeros(0)
    Aswitch = kron(ymarkov_combined, speye(Float64, I))

    profshare, r, lumptransfer = calculate_eqcond_equil_objects(vars_SS, params, zz, w, TFP, MP)

    ## HJB
    h0 = copy(h_SS)
    h  = copy(h_SS)

    util, income, labor = construct_household_problem_functions(V, w, params)

    Vaf, Vab, cf, hf, cb, hb = construct_initial_diff_matrices(V, Vaf, Vab, income, labor, h,
                                                               h0, zz, profshare,
                                                               lumptransfer, amax, amin,
                                                               coefrra, r, dazf, dazb,
                                                               maxhours)

    cf, hf, cb, hb, c0, h0 = hours_iteration(income, labor, zz, profshare, lumptransfer, aa,
                                             coefrra, r, cf, hf, cb, hb, c0, h0, maxhours, niter_hours)

    for i in eachindex(h0)
        c0[i] = income(h0[i], zz[i], profshare[i], lumptransfer, r, aa[i])

        sf[i] = income(hf[i], zz[i], profshare[i], lumptransfer, r, aa[i]) - cf[i]
        sb[i] = income(hb[i], zz[i], profshare[i], lumptransfer, r, aa[i]) - cb[i]
    end

    Vaf[I, :] = cf[I, :] .^ (-coefrra)
    Vab[1, :] = cb[1, :] .^ (-coefrra)
    Va0 = (c0).^(-coefrra)

    V, A, u, h, c, s  = upwind(V, util, Aswitch, cf, cb, c0, hf, hb, h0, sf, sb, Vaf, Vab, Va0,
                               dazf, dazb)

    ## Collect dynamics
    v_residual = calculate_residuals(vars_SS, params, u, A, V, VDot, VEErrors, w, TFP,
                                     inflationDot, inflationEError, azdelta_mat, g, MP,
                                     MPDot, MPShock, aa, zz, azdelta, r, gDot, s, h, c)

    #    return eqm conditions with v_residual
    return v_residual
end

# Solve out the static constraints of the model
function solve_static_constraints(Γ0::Matrix, Γ1::Matrix, Ψ::Matrix, Π::Matrix, C::Matrix)
    # Identify zero rows of Γ0
    n = size(Γ0, 1)
    redundant = find(sum(abs.(Γ0), 2) .== 0)
    keep = setdiff(collect(1:n), redundant)
    n_keep = length(keep)

    redundant_states_inv = zeros(n, n_keep)
    redundant_states_inv[keep, :] = eye(n_keep)
    redundant_states_inv[redundant, :] = -Γ1[redundant, redundant]\Γ1[redundant, keep]

    redundant_states = zeros(n_keep, n)
    redundant_states[:, keep] = eye(n_keep)

    Γ0 = redundant_states * Γ0 * redundant_states_inv
    Γ1 = redundant_states * Γ1 * redundant_states_inv
    Γ1 = Γ0\Γ1
    Γ0_red = Γ0\redundant_states
    Ψ  = Γ0_red * Ψ
    Π  = Γ0_red * Π
    C  = Γ0_red * C

    return Γ0, Γ1, Ψ, Π, C, redundant_states_inv
end