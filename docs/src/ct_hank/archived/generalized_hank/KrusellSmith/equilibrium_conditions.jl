# Set up model in canonical form
function equilibrium_conditions(varsSS::OrderedDict{Symbol, Any},
                       grids::Dict{Symbol, Any},
                       params::Dict{Symbol, Float64},
                       approx_params::Dict{Symbol, Any})

    nVars = grids[:nVars]
    nEErrors = grids[:nEErrors]
    x = zeros(Float64, 2*nVars + nEErrors + 1)
    function f{T<:Real}(x::Vector{T}; vars_ss = deepcopy(varsSS))
        v_residual = get_residuals(x, vars_ss, grids, params, approx_params)
        return vcat(values(v_residual)...)
    end

    derivs = ForwardDiff.jacobian(f, x)
    Γ1 = -derivs[:,1:nVars]
    Γ0 = derivs[:,nVars+1:2*nVars]
    Π = -derivs[:,2*nVars+1:2*nVars+nEErrors]
    Ψ = -derivs[:,2*nVars+nEErrors+1:2*nVars+ nEErrors + 1]
    C = zeros(nVars, 1)

    canonical_form = Dict{Symbol, Any}()
    canonical_form[:Γ1] = Γ1
    canonical_form[:Γ0] = Γ0
    canonical_form[:Π] = Π
    canonical_form[:Ψ] = Ψ
    canonical_form[:C] = C

    return canonical_form
end


function get_residuals(varsSS::OrderedDict{Symbol, Any},
                       grids::Dict{Symbol, Any},
                       params::Dict{Symbol, Float64},
                       approx_params::Dict{Symbol, Any})
    nVars = grids[:nVars]
    nErrors = grids[:nEErrors]
    x = zeros(2*nVars + nEErrors + 1)
    return get_residuals(x, varsSS, grids, params, approx_params)
end

function get_residuals{T<:Real}(x::Vector{T},
                       varsSS::OrderedDict{Symbol, Any},
                       grids::Dict{Symbol, Any},
                       params::Dict{Symbol, Float64},
                       approx_params::Dict{Symbol, Any})

    # Read in steady state vars, params, and grids
    γ = params[:γ]
    ρ = params[:ρ]
    δ = params[:δ]
    α = params[:α]
    σ_tfp = params[:σ_tfp]
    ρ_tfp = params[:ρ_tfp]
    μ = params[:μ]
    τ = params[:τ]
    I = grids[:I]
    amin = grids[:amin]
    amax = grids[:amax]
    a = grids[:a]
    da = grids[:da]
    aa = grids[:aa]
    zz = grids[:zz]
    Aswitch = grids[:Aswitch]
    aaa = grids[:aaa]
    zzz = grids[:zzz]
    zAvg = grids[:zAvg]
    nVars = grids[:nVars]
    nEErrors = grids[:nEErrors]

    # Unpack vars
    V = x[1:2*I] + varsSS[:VSS]
    g = x[2*I+1:4*I-1] + varsSS[:ggSS]
    g_end=1/da-sum(g) # ensures that distribution integrates to 1 when multiplied by da
    logAggregateTFP = x[4*I]
    KHat = x[4*I+1] + varsSS[:KSS]
    rHat = x[4*I+2] + varsSS[:rSS]
    wHat = x[4*I+3] + varsSS[:wSS]
    output = x[4*I+4] + varsSS[:YSS]
    C = x[4*I+5] + varsSS[:CSS]
    investment = x[4*I+6] + varsSS[:ISS]
    IfSS = varsSS[:IfSS]
    IbSS = varsSS[:IbSS]
    I0SS = varsSS[:I0SS]

    V = reshape(V,I,2)
    VDot = x[nVars+1:nVars+2*I]
    gDot = x[nVars+2*I+1:nVars+4*I-1]
    logAggregateTFPDot = x[nVars+4*I]
    VEErrors = x[2*nVars+1:2*nVars+2*I]
    aggregateTFPShock = x[2*nVars+nEErrors+1]

    # Initialize other variables, using vars to ensure everything is a dual number
    dVf = copy(V) #what if I make this a copy instead?
    dVb = copy(V)

    K = sum(aaa .* [g;g_end] * da)
    r = exp(logAggregateTFP) * α * (KHat ^ (α - 1)) * (zAvg ^ (1 - α)) - δ
    w = exp(logAggregateTFP) * (1 - α) * (KHat ^ α) * (zAvg ^ (-α))

    #----------------------------------------------------------------
    # Compute one iteration of HJB Equation
    #----------------------------------------------------------------
    c0 = w * ((1 - τ) * zz + μ * (1 - zz)) + r * aa

    # Compute forward difference
    dVf[1:I-1,:] = (V[2:I,:]-V[1:I-1,:])/da
    dVf[I,:] = c0[I,:] .^ (-γ) #will never be used, but impose state constraint a<=amax just in case

    # Compute backward difference
    dVb[2:I,:] = (V[2:I,:]-V[1:I-1,:])/da
    dVb[1,:] = c0[1,:] .^ (-γ) #state constraint boundary condition

    # Compute consumption and savings with forward difference
    cf = dVf.^(-1/γ)
    ssf = c0 - cf

    # Compute consumption and savings with backward difference
    cb = dVb.^(-1/γ)
    ssb = c0 - cb

    # Compute consumption and derivative of value function for no drift
    dV0 = c0.^(-γ)

    # Compute upwind difference
    dV_Upwind = dVf.*IfSS + dVb.*IbSS + dV0.*I0SS
    c = dV_Upwind.^(-1/γ)
    u = c.^(1-γ)./(1-γ)
    savings = c0 - c

    # Construct A matrix
    X = -ssb.*IbSS/da
    Y = -ssf.*IfSS/da + ssb.*IbSS/da
    Z = ssf.*IfSS/da

    X[1,:]=0
    lowdiag=reshape(X,2*I,1)
    Z[I,:]=0
    A = spdiagm(tuple(reshape(Y,2*I,1)),0,2*I,2*I) + spdiagm(lowdiag[2:2*I],-1,2*I,2*I) + spdiagm(tuple(reshape(Z,1,2*I)[1:(2*I-1)]),1,2*I,2*I) + Aswitch

    #----------------------------------------------------------------
    # Compute equilibrium conditions
    #----------------------------------------------------------------
    # HJB Equation
    hjbResidual = reshape(u,2*I,1) + A * reshape(V,2*I,1) + VDot + VEErrors - ρ * reshape(V,2*I,1)

     # KFE
    gIntermediate = A' * [g;g_end]
    gResidual = gDot - gIntermediate[1:2*I-1,1]

    # Aggregates
    kResidual = K - KHat
    rResidual = r - rHat
    wResidual = w - wHat
    yResidual = output - exp(logAggregateTFP) * (sum(aaa .* [g;g_end] * da) ^ α) * (zAvg ^ (1 - α))
    cResidual = C - sum(c[:] .* [g;g_end] * da)
    iResidual = investment - sum((savings[:] + δ * aaa) .* [g;g_end] * da)

    # Law of motion for aggregate shocks
    tfpResidual = logAggregateTFPDot + (1 - ρ_tfp) * logAggregateTFP - σ_tfp * aggregateTFPShock

    # Store residuals
    v_residual = OrderedDict{Symbol, Any}()
    v_residual[:hjbResidual] = hjbResidual
    v_residual[:gResidual] = gResidual
    v_residual[:tfpResidual] = tfpResidual
    v_residual[:kResidual] = kResidual
    v_residual[:rResidual] = rResidual
    v_residual[:wResidual] = wResidual
    v_residual[:yResidual] = yResidual
    v_residual[:cResidual] = cResidual
    v_residual[:iResidual] = iResidual

    return v_residual
end
