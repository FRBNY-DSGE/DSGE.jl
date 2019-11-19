# Computes (non-stochastic) steady-state values in the two-asset Krusell Smith HANK model


function compute_steady_state(grids::Dict{Symbol, Any}, params::Dict{Symbol, Float64},
                              init_params::Dict{Symbol, Float64}, approx_params::Dict{Symbol, Any})

    # Read in parameters
    γ = params[:γ]
    ρ = params[:ρ]
    δ = params[:δ]
    α = params[:α]
    σ_tfp = params[:σ_tfp]
    ρ_tfp = params[:ρ_tfp]
    μ = params[:μ] # UI replacement rate
    τ = params[:τ] # labor income tax

    # Read in grids
    lla = grids[:lla] # transition probabilities in idiosyncratic income shocks
    z = grids[:z]
    a = grids[:a]
    da = grids[:da]
    dz = grids[:dz]
    amax = grids[:amax]
    amin = grids[:amin]
    zAvg = grids[:zAvg] # average labor supply
    I = grids[:I]
    aa = grids[:aa]
    zz = grids[:zz]
    Aswitch = grids[:Aswitch]

    # Read in initial rates
    r0 = init_params[:r0]
    rmin = init_params[:rmin]
    rmax = init_params[:rmax]

    # Read in  approximation parameters
    Ir = approx_params[:Ir]
    crit = approx_params[:crit]
    Δ = approx_params[:Δ]
    crit_S = approx_params[:crit_S]
    maxit_HJB = approx_params[:maxit_HJB]

    # Initialize steady state variables
    vars_SS = OrderedDict{Symbol, Any}()

    # Initialize guesses and difference matrices
    dVf = zeros(I, 2)
    dVb = zeros(I, 2)
    dV_upwind = zeros(I, 2)
    c = zeros(I, 2)
    u = zeros(I, 2)
    g = zeros(I, 2)
    r = r0
    KD = (((α) / (r + δ)) ^ (1 / (1 - α))) * zAvg
    w = (1 - α) * (KD ^ α) * ((zAvg) ^ (-α))
    v0 = zeros(I, 2)
    v0[:, 1] = (w * μ * (1 - z[1]) + r.*a).^(1 - γ)/(1 - γ)/ρ
    v0[:, 2] = (w * (1 - τ) * z[2] + r.*a).^(1 - γ)/(1 - γ)/ρ

    If = zeros(I, 2)
    Ib = zeros(I, 2)
    I0 = zeros(I,2)
    A = spzeros(2*I, 2*I)

    # Iterate to find steady state interest rate
    for ir = 1:Ir
        # Get capital demand and wages
        KD = (((α) / (r + δ)) ^ (1 / (1 - α))) * zAvg
        w = (1 - α) * (KD ^ α) * ((zAvg) ^ (-α))

        v = v0

        # Solve for value function given interest rate
        for n = 1:maxit_HJB
            # Grab guess
            V = v

            # Compute forward difference
            dVf[1:I - 1, :] = (V[2:I, :] - V[1:(I - 1), :])/da
            dVf[I, :] = (w * ((1 - τ) * z + μ * (1 - z)) + r.*amax).^(-γ)

            # Compute backward difference
            dVb[2:I, :] = (V[2:I, :] - V[1:(I - 1), :])/da
            dVb[1, :] = (w * ((1 - τ) * z + μ * (1 - z)) + r.*amin).^(-γ)

            # Compute consumption and savings with forward difference
            cf = dVf.^(-1/γ)
            ssf = w*((1 - τ) * zz + μ * (1 - zz)) + r.*aa - cf

            # Compute consumption and savings with backward difference
            cb = dVb.^(-1/γ)
            ssb = w*((1 - τ) * zz + μ * (1 - zz)) + r.*aa - cb

            # Compute consumption and derivative of value function for no drift
            c0 = w*((1 - τ) * zz + μ * (1 - zz)) + r.*aa
            dV0 = c0.^(-γ)

            # Compute upwind on value function
            dV_upwind, If, Ib, I0 = upwind_value_function(dVf, dVb, dV0, ssf, ssb)
            c = dV_upwind.^(-1 / γ)
            u = c.^(1 - γ)/(1 - γ)
            savingsSS = w * ((1 - τ) * zz + μ * (1 - zz)) + r.*aa - c

            # Compute A matrix
            A = upwind_matrix(Aswitch, ssf, ssb, da, I, 2)

            # Solve for value function using upwind scheme
            V_stacked = solve_hjb(A, ρ, Δ, u, V)
            V = reshape(V_stacked, I, 2)

            # Update values and check convergence
            Vchange = V - v
            v = V
            if maximum(abs.(Vchange)) < crit
                break
            end
        end
        v0 = v

        # Solve for stationary distribution
        g = solve_KFE(A, grids)

        # Compute objects from this iteration
        KS = g[:,1]' * a * da + g[:,2]' * a * da
        S = KS - KD

        # Update interest rate
        if S > crit_S
            # Excess Supply
            rmax = r
            r = 0.5 * (r + rmin)
        elseif S < -crit_S
            # Excess Demand
            rmin = r
            r = 0.5 * (r + rmax)
        elseif abs(S) < crit_S
            K = KS
            # grab SS values
            rSS = r
            wSS = w
            KSS = K
            ASS = A
            uSS = u
            cSS = c
            VSS = v
            gSS = g
            dVUSS = dV_upwind
            dVfSS = dVf
            dVbSS = dVb
            IfSS = If
            IbSS = Ib
            I0SS = I0

            vars_SS[:VSS] = reshape(VSS, 2*I, 1)
            ggSS = reshape(gSS, 2*I, 1)
            vars_SS[:ggSS] = ggSS[1:2*I-1]
            vars_SS[:logAggregateTFP] = 0
            vars_SS[:KSS] = KSS
            vars_SS[:rSS] = rSS
            vars_SS[:wSS] = wSS
            vars_SS[:YSS] = (KSS ^ α) * (zAvg ^ (1 - α))
            CSS = sum(cSS[:] .* gSS[:] * da)
            vars_SS[:CSS] = CSS
            vars_SS[:ISS] = δ * KSS
            vars_SS[:IfSS] = IfSS
            vars_SS[:IbSS] = IbSS
            vars_SS[:I0SS] = I0SS

            break
        end
    end

   return vars_SS
end

