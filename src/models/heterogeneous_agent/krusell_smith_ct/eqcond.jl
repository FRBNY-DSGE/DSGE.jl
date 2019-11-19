"""
```
eqcond(m::KrusellSmithCT)
```

Expresses the equilibrium conditions in canonical form using Γ0, Γ1, C, Ψ, and Π matrices.

### Outputs

* `Γ0` (`n_vars` x `n_vars`) holds coefficients of current time states.
* `Γ1` (`n_vars` x `n_vars`) holds coefficients of lagged states.
* `C`  (`n_vars` x `1`) is a vector of constants
* `Ψ`  (`n_vars` x `n_shocks_exogenous`) holds coefficients of iid shocks.
* `Π`  (`n_vars` x `n_states_expectational`) holds coefficients of expectational states.
"""
function eqcond(m::KrusellSmithCT)

    # Read in steady state vars, params, and grids
    nstates = n_states(m)
    x       = zeros(2 * nstates + n_shocks_expectational(m) + n_shocks_exogenous(m))

    γ::Float64     = m[:γ].value
    ρ::Float64     = m[:ρ].value
    δ::Float64     = m[:δ].value
    α::Float64     = m[:α].value

    σ_tfp::Float64 = m[:σ_tfp].value
    ρ_tfp::Float64 = m[:ρ_tfp].value

    λ1::Float64    = m[:λ1].value
    λ2::Float64    = m[:λ2].value

    μ::Float64     = m[:μ].value
    τ::Float64     = m[:τ].value
    I::Int64       = get_setting(m, :I)
    J::Int64       = get_setting(m, :J)
    IJ             = I*J

    amin::Float64  = get_setting(m, :amin)
    amax::Float64  = get_setting(m, :amax)
    da::Float64    = get_setting(m, :da)

    If_ss::Array{Float64,2} = m[:If_ss].value
    Ib_ss::Array{Float64,2} = m[:Ib_ss].value
    I0_ss::Array{Float64,2} = m[:I0_ss].value

    a::Array{Float64,1}     = get_setting(m, :a)
    aa::Array{Float64,2}    = get_setting(m, :aa)
    aaa::Array{Float64,2}   = get_setting(m, :aaa)
    z::Array{Float64,1}     = get_setting(m, :z)
    zz::Array{Float64,2}    = get_setting(m, :zz)
    zzz::Array{Float64,2}   = get_setting(m, :zzz)

    z_avg::Float64          = (λ1 * z[2] + λ2 * z[1]) / (λ1 + λ2)

    V_ss  = m[:V_ss].value
    gg_ss = m[:gg_ss].value
    K_ss  = m[:K_ss].value
    r_ss  = m[:r_ss].value
    w_ss  = m[:w_ss].value
    Y_ss  = m[:Y_ss].value
    C_ss  = m[:C_ss].value
    I_ss  = m[:I_ss].value

    A_switch = [-SparseMatrixCSC{Float64}(LinearAlgebra.I, I, I) * λ1 SparseMatrixCSC{Float64}(LinearAlgebra.I, I, I) * λ1; SparseMatrixCSC{Float64}(LinearAlgebra.I, I, I) * λ2 -SparseMatrixCSC{Float64}(LinearAlgebra.I, I, I) * λ2]

    function get_residuals(x::Vector{T}) where {T<:Real}

        # Unpack vars and convert all values to dual types (see JuliaDiff)
        V_vec             = x[1:2*I]       + V_ss
        gg                = x[2*I+1:4*I-1] + gg_ss
        # Ensures that distribution integrates to 1 when multiplied by da
        g_end             = 1/da - sum(gg)
        g                 = [gg;g_end]
        g_dot_da          = g * da
        log_aggregate_tfp = x[4*I]
        K_hat             = x[4*I+1] + K_ss
        r_hat             = x[4*I+2] + r_ss
        w_hat             = x[4*I+3] + w_ss
        output            = x[4*I+4] + Y_ss
        C_cur             = x[4*I+5] + C_ss
        investment        = x[4*I+6] + I_ss

        # Set up system of time derivatives
        V                     = reshape(V_vec, I, J)
        V_dot                 = x[  nstates + 1 : nstates + 2*I]
        g_dot                 = x[  nstates + 2*I + 1 : nstates + 4*I - 1]
        log_aggregate_tfp_dot = x[  nstates + 4*I]
        V_exp_errors          = x[2*nstates + 1 : 2*nstates + 2*I] # expected errors in the value function
        aggregate_tfp_shock   = x[end]

        # Instantiate variables for future use
        c = similar(V_vec)
        u = similar(V_vec)
        X = similar(V)
        Y = similar(V)
        Z = similar(V)

        # Get steady state values of K, r, w incorporating steady state aggregate_tfp level
        K = sum(aaa .* g_dot_da)
        r = exp(log_aggregate_tfp) * α * (K_hat ^ (α - 1)) * (z_avg ^ (1 - α)) - δ
        w = exp(log_aggregate_tfp) * (1 - α) * (K_hat ^ α) * (z_avg ^ (-α))

        #----------------------------------------------------------------
        # Compute one iteration of HJB Equation
        #----------------------------------------------------------------
        c0 = w * ((1 - τ) * zz + μ * (1 .- zz)) + r * aa

        for j=1:J, i=1:I
            # Compute consumption and derivative of value function for no drift
            dV0 = c0[i,j] ^ (-γ)

            if i==I
                dVf = dV0 # Will never be used, but impose state constraint a<=amax just in case
                dVb = (V[i,j] - V[i-1,j]) / da
            elseif i==1
                dVb = dV0                      # State constraint boundary condition
                dVf = (V[i+1,j] - V[i,j]) / da
            else
                dVf = (V[i+1,j] - V[i,j]) / da # Compute forward difference
                dVb = (V[i,j] - V[i-1,j]) / da # Compute backward difference
            end

            # Compute consumption and savings with forward difference
            ssf = c0[i,j] - dVf ^ (-1.0 / γ)

            # Compute consumption and savings with backward difference
            ssb = c0[i,j] - dVb ^ (-1.0 / γ)

            # Compute consumption and derivative of value function for no drift, then compute upwind difference
            c[(j-1)*I + i] = (dVf * If_ss[i,j] + dVb * Ib_ss[i,j] + dV0 * I0_ss[i,j]) ^ (-1.0 / γ)
            u[(j-1)*I + i] = (c[(j-1)*I + i] ^ (1.0 - γ)) / (1.0 - γ)

            # Construct A matrix
            X[i,j] = -ssb * (Ib_ss[i,j] / da)
            Z[i,j] =  ssf * (If_ss[i,j] / da)
            Y[i,j] = -X[i,j] - Z[i,j]

            # Note: Need to handle complex 0 case; should check if it ever gets there
            if i==I
                Z[i,j] = 0.
            elseif i==1
                X[i,j] = 0.
            end
        end

        A = spdiagm(0 => reshape(Y,IJ), -1 => reshape(X,IJ)[2:IJ], 1 => reshape(Z,IJ)[1:IJ-1]) + A_switch

        #----------------------------------------------------------------
        # Compute residuals of equilibrium conditions
        #----------------------------------------------------------------

        # HJB Equation
        hjb_residual = u + A * V_vec + V_dot + V_exp_errors - ρ * V_vec

        # KFE
        g_residual   = g_dot - (A' * g)[1:IJ-1, 1]

        # Law of motion for aggregate shocks
        tfp_residual = log_aggregate_tfp_dot + (1 - ρ_tfp) * log_aggregate_tfp - σ_tfp * aggregate_tfp_shock

        # Aggregates
        k_residual = K - K_hat
        r_residual = r - r_hat
        w_residual = w - w_hat
        y_residual = output - exp(log_aggregate_tfp) * (K ^ α) * (z_avg ^ (1 - α))
        c_residual = C_cur - sum(c .* g_dot_da)
        i_residual = investment - sum((reshape(c0, IJ) - c + δ * aaa) .* g_dot_da)

        return  [hjb_residual; g_residual; tfp_residual; k_residual; r_residual; w_residual;
                 y_residual; c_residual; i_residual]
    end

    #----------------------------------------------------------------
    # Calculate Jacobian, unpack output
    #----------------------------------------------------------------
    derivs = ForwardDiff.jacobian(get_residuals, x)

    Γ1::Array{Float64,2} = -derivs[:, 1 : nstates]
    Γ0::Array{Float64,2} =  derivs[:, nstates + 1 : 2 * nstates]
    Π::Array{Float64,2}  = -derivs[:, 2 * nstates + 1 : 2 * nstates + n_shocks_expectational(m)]
    Ψ::Array{Float64,2}  = -derivs[:, 2 * nstates + n_shocks_expectational(m) + 1 : end]
    C::Array{Float64,1}  = zeros(nstates)

    return Γ0, Γ1, Ψ, Π, C
end