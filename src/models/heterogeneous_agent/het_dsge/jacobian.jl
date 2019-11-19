function jacobian(m::HetDSGE)
    truncate_distribution!(m)

    # Load in endogenous state and eq cond indices
    endo = augment_model_states(m.endogenous_states_unnormalized,
                         n_model_states_unnormalized(m))
    eq   = m.equilibrium_conditions

    # Load in parameters, steady-state parameters, and grids
    r::Float64     = m[:r].value
    α::Float64     = m[:α].value
    H::Float64     = m[:H].value
    δ::Float64     = m[:δ].value
    μ_sp::Float64  = m[:μ_sp].value
    ρ_sp::Float64  = m[:ρ_sp].value
    σ_sp::Float64  = m[:σ_sp].value
    γ::Float64     = m[:γ].value
    g::Float64     = m[:g].value
    η::Float64     = m[:η].value
    ρ_B::Float64    = m[:ρ_B].value
    ρ_G::Float64    = m[:ρ_G].value
    ρ_z::Float64    = m[:ρ_z].value
    ρ_μ::Float64    = m[:ρ_μ].value
    ρ_lamw::Float64 = m[:ρ_lamw].value
    ρ_lamf::Float64 = m[:ρ_lamf].value
    ρ_mon::Float64  = m[:ρ_mon].value
    spp::Float64   = m[:spp].value
    ϕh::Float64    = m[:ϕh].value
    ρ_R::Float64    = m[:ρR].value
    ψπ::Float64    = m[:ψπ].value
    ψy::Float64    = m[:ψy].value
    κ_p::Float64  = m[:κ_p].value
    κ_w::Float64  = m[:κ_w].value


    R = 1 + r

    ell::Vector{Float64}  = m[:lstar].value
    c::Vector{Float64}    = m[:cstar].value
    μ::Vector{Float64}    = m[:μstar].value
    β::Float64            = m[:βstar].value

    T::Float64     = m[:Tstar].value
    ω::Float64     = m[:ωstar].value
    xstar::Float64 = m[:xstar].value
    ystar::Float64 = m[:ystar].value
    Rk::Float64    = m[:Rkstar].value
    kstar::Float64 = m[:kstar].value


    xgrid::Vector{Float64} = m.grids[:xgrid].points
    xwts::Vector{Float64}  = m.grids[:xgrid].weights
    sgrid::Vector{Float64} = m.grids[:sgrid].points
    swts::Vector{Float64}  = m.grids[:sgrid].weights
    fgrid::Matrix{Float64} = m.grids[:fgrid]
    xswts = kron(swts,xwts)

    zlo::Float64 = get_setting(m, :zlo)
    zhi::Float64 = get_setting(m, :zhi)

    nx::Int = get_setting(m, :nx)
    ns::Int = get_setting(m, :ns)
    nxns = nx*ns

    qp(z) = dmollifier_hetdsge(z, zhi, zlo)
    qfunction_hetdsge(x) = mollifier_hetdsge(x, zhi, zlo) #/sumz

    unc = 1 ./ ell .<= repeat(xgrid,ns) .+ η

    dF1_dELL, dF1_dRZ, dF1_dELLP, dF1_dWHP, dF1_dTTP, ee =
        euler_equation_hetdsge_lag(nx, ns, qp, qfunction_hetdsge, xgrid, sgrid, fgrid, unc, xswts,
                               R, γ, β, η, ell, T, ω, H)

    # KF Equation
    dF2_dWH, dF2_dRZ, dF2_dTT,dF2_dELL, bigΨ, dF2_dM =
        kolmogorov_fwd_hetdsge(nx, ns, qfunction_hetdsge, qp, xgrid, sgrid, fgrid, unc, xswts,
                               R, γ, ell, μ, η, T, ω, H, ee)

    # Market clearing, lambda function
    c = min.(1 ./ ell,repeat(xgrid,ns).+η)
    lam = (xswts.*μ)'*(1 ./ c) # average marginal utility which the union uses to set wages
    ϕ = lam*ω/(H^ϕh) # now that we know lam in steady state, choose disutility to target hours H

    setup_indices!(m)
    normalize_model_state_indices!(m)

    nvars = get_setting(m, :nvars)
    # Make the Jacobian
    JJ = zeros(nvars, 2*nvars)

    # Euler equation
    JJ[eq[:eq_euler],endo[:l′_t]] = dF1_dELLP
    JJ[eq[:eq_euler],endo[:z′_t]]   = -dF1_dRZ
    JJ[eq[:eq_euler],endo[:w′_t]]   = dF1_dWHP
    JJ[eq[:eq_euler],endo[:L′_t]]  = dF1_dWHP
    JJ[eq[:eq_euler],endo[:t′_t]]  = dF1_dTTP
    JJ[eq[:eq_euler],endo[:b_t]]    = ell
    JJ[eq[:eq_euler],endo[:l_t]]  = dF1_dELL
    JJ[eq[:eq_euler],endo[:R_t]]   = dF1_dRZ

    # KF eqn
    JJ[eq[:eq_kolmogorov_fwd],endo[:kf′_t]]   = -Matrix{Float64}(I, nxns, nxns)
    JJ[eq[:eq_kolmogorov_fwd],endo[:kf_t]]   = dF2_dM
    JJ[eq[:eq_kolmogorov_fwd],endo[:l_t]] = dF2_dELL
    JJ[eq[:eq_kolmogorov_fwd],endo[:R_t]] = dF2_dM*dF2_dRZ
    JJ[eq[:eq_kolmogorov_fwd],endo[:z_t]]    = -dF2_dM*dF2_dRZ
    JJ[eq[:eq_kolmogorov_fwd],endo[:w_t]]    = dF2_dM*dF2_dWH
    JJ[eq[:eq_kolmogorov_fwd],endo[:L_t]]   = dF2_dM*dF2_dWH
    JJ[eq[:eq_kolmogorov_fwd],endo[:t_t]]   = dF2_dM*dF2_dTT

    # mkt clearing
    JJ[first(eq[:eq_market_clearing]),first(endo[:y_t])]   = ystar/g
    JJ[first(eq[:eq_market_clearing]),first(endo[:g_t])]   = -ystar/g
    JJ[first(eq[:eq_market_clearing]),first(endo[:I_t])]   = -xstar
    JJ[first(eq[:eq_market_clearing]), endo[:l_t]] = (μ .*unc.*xswts.*c)'
    JJ[first(eq[:eq_market_clearing]), endo[:kf_t]] = -(xswts.*c)' # note, now we linearize
    JJ[first(eq[:eq_market_clearing]),first(endo[:R_t])]   = -(xswts.*c)'*dF2_dRZ
    JJ[first(eq[:eq_market_clearing]),first(endo[:z_t])]   = (xswts.*c)'*dF2_dRZ
    JJ[first(eq[:eq_market_clearing]),first(endo[:w_t])]   = -(xswts.*c)'*dF2_dWH
    JJ[first(eq[:eq_market_clearing]),first(endo[:L_t])]   = -(xswts.*c)'*dF2_dWH
    JJ[first(eq[:eq_market_clearing]),first(endo[:t_t])]   = -(xswts.*c)'*dF2_dTT

    # lambda = average marginal utility
    JJ[first(eq[:eq_lambda]),first(endo[:margutil_t])] = lam
    JJ[first(eq[:eq_lambda]),endo[:kf_t]]   = -(xswts./c)' # note, now we linearize
    JJ[first(eq[:eq_lambda]),first(endo[:R_t])]   = -(xswts./c)'*dF2_dRZ
    JJ[first(eq[:eq_lambda]),first(endo[:z_t])]   = (xswts./c)'*dF2_dRZ
    JJ[first(eq[:eq_lambda]),first(endo[:w_t])]   = -(xswts./c)'*dF2_dWH
    JJ[first(eq[:eq_lambda]),first(endo[:L_t])]   = -(xswts./c)'*dF2_dWH
    JJ[first(eq[:eq_lambda]),first(endo[:t_t])]   = -(xswts./c)'*dF2_dTT
    JJ[first(eq[:eq_lambda]),endo[:l_t]] = -(xswts.*unc.*μ./c)'

    # transfer
    JJ[first(eq[:eq_transfers]),first(endo[:t_t])]  = T
    JJ[first(eq[:eq_transfers]),first(endo[:capreturn_t])]  = -Rk*kstar
    JJ[first(eq[:eq_transfers]),first(endo[:k_t])]  = -Rk*kstar
    JJ[first(eq[:eq_transfers]),first(endo[:z_t])]  = Rk*kstar
    JJ[first(eq[:eq_transfers]),first(endo[:I_t])]   = xstar
    JJ[first(eq[:eq_transfers]),first(endo[:mc_t])]  = ystar
    JJ[first(eq[:eq_transfers]),first(endo[:y_t])]   = (1-1/g)*ystar
    JJ[first(eq[:eq_transfers]),first(endo[:g_t])]   = (ystar/g)

    # investment
    JJ[first(eq[:eq_investment]),first(endo[:Q_t])]  = 1.
    JJ[first(eq[:eq_investment]),first(endo[:μ_t])] = 1.
    JJ[first(eq[:eq_investment]),first(endo[:I′_t])] = spp*(exp(3*γ))/R
    JJ[first(eq[:eq_investment]),first(endo[:z′_t])] = spp*(exp(3*γ))/R
    JJ[first(eq[:eq_investment]),first(endo[:I_t])]  = -spp*(exp(3*γ))/R - spp*exp(2*γ)
    JJ[first(eq[:eq_investment]),first(endo[:I_t1])] = spp*exp(2*γ)
    JJ[first(eq[:eq_investment]),first(endo[:z_t])]  = -spp*exp(2*γ)

    # tobin's q
    JJ[first(eq[:eq_tobin_q]),first(endo[:R_t])]  = R
    JJ[first(eq[:eq_tobin_q]),first(endo[:Q_t])]   = R
    JJ[first(eq[:eq_tobin_q]),first(endo[:capreturn′_t])] = -Rk
    JJ[first(eq[:eq_tobin_q]),first(endo[:Q′_t])]  = -(1-δ)

    # capital accumulation
    JJ[first(eq[:eq_capital_accumulation]),first(endo[:k′_t])] = 1.
    JJ[first(eq[:eq_capital_accumulation]),first(endo[:k_t])]  = -(1-δ)
    JJ[first(eq[:eq_capital_accumulation]),first(endo[:z_t])]   = (1-δ)
    JJ[first(eq[:eq_capital_accumulation]),first(endo[:μ_t])]  = -xstar/kstar
    JJ[first(eq[:eq_capital_accumulation]),first(endo[:I_t])]   = -xstar/kstar

    # wage phillips curve
    JJ[first(eq[:eq_wage_phillips]),first(endo[:π_w_t])]  = -1.
    JJ[first(eq[:eq_wage_phillips]),first(endo[:λ_w_t])] = 1. #(ϕ*H^ϕh)/Φw
    JJ[first(eq[:eq_wage_phillips]),first(endo[:L_t])]   = κ_w*ϕh #(ϕ*(H^ϕh)*(1+lamw)/lamw*Φw)*ϕh
    JJ[first(eq[:eq_wage_phillips]),first(endo[:margutil_t])]  = -κ_w #-(ϕ*(H^ϕh)*(1+lamw)/lamw*Φw)
    JJ[first(eq[:eq_wage_phillips]),first(endo[:w_t])]    = -κ_w #-(ϕ*(H^ϕh)*(1+lamw)/lamw*Φw)
    JJ[first(eq[:eq_wage_phillips]),first(endo[:π_w′_t])]  = β

    # price phillips curve
    JJ[first(eq[:eq_price_phillips]),first(endo[:π_t])]   = -1.
    JJ[first(eq[:eq_price_phillips]),first(endo[:mc_t])]   = κ_p #(1+lamf)/(lamf*Φp)
    JJ[first(eq[:eq_price_phillips]),first(endo[:λ_f_t])] = 1. #1/Φp
    JJ[first(eq[:eq_price_phillips]),first(endo[:π′_t])]  = 1/R

    # marginal cost
    JJ[first(eq[:eq_marginal_cost]),first(endo[:mc_t])] = 1.
    JJ[first(eq[:eq_marginal_cost]),first(endo[:w_t])]  = -(1-α)
    JJ[first(eq[:eq_marginal_cost]),first(endo[:capreturn_t])] = -α

    # gdp
    JJ[first(eq[:eq_gdp]),first(endo[:y_t])]  = 1.
    JJ[first(eq[:eq_gdp]),first(endo[:z_t])]  = α
    JJ[first(eq[:eq_gdp]),first(endo[:k_t])] = -α
    JJ[first(eq[:eq_gdp]),first(endo[:L_t])] = -(1-α)

    # optimal k/l ratio
    JJ[first(eq[:eq_optimal_kl]),first(endo[:capreturn_t])] = 1.
    JJ[first(eq[:eq_optimal_kl]),first(endo[:w_t])]  = -1.
    JJ[first(eq[:eq_optimal_kl]),first(endo[:L_t])] = -1.
    JJ[first(eq[:eq_optimal_kl]),first(endo[:k_t])] = 1.
    JJ[first(eq[:eq_optimal_kl]),first(endo[:z_t])]  = -1.

    # taylor rule
    JJ[first(eq[:eq_taylor]),first(endo[:i_t])]   = -1.
    JJ[first(eq[:eq_taylor]),first(endo[:i_t1])]  = ρ_R
    JJ[first(eq[:eq_taylor]),first(endo[:π_t])]   = (1-ρ_R)*ψπ
    JJ[first(eq[:eq_taylor]),first(endo[:y_t])]    = (1-ρ_R)*ψy
    JJ[first(eq[:eq_taylor]),first(endo[:y_t1])]  = -(1-ρ_R)*ψy
    JJ[first(eq[:eq_taylor]),first(endo[:z_t])]    = (1-ρ_R)*ψy
    JJ[first(eq[:eq_taylor]),first(endo[:rm_t])] = 1.

    # fisher eqn
    JJ[first(eq[:eq_fisher]),first(endo[:R_t])]  = 1.
    JJ[first(eq[:eq_fisher]),first(endo[:π′_t])] = 1.
    JJ[first(eq[:eq_fisher]),first(endo[:i_t])]  = -1.

    # wage inflation
    JJ[first(eq[:eq_nominal_wage_inflation]),first(endo[:π_w_t])] = 1.
    JJ[first(eq[:eq_nominal_wage_inflation]),first(endo[:π_t])]  = -1.
    JJ[first(eq[:eq_nominal_wage_inflation]),first(endo[:z_t])]   = -1.
    JJ[first(eq[:eq_nominal_wage_inflation]),first(endo[:w_t])]   = -1.
    JJ[first(eq[:eq_nominal_wage_inflation]),first(endo[:w_t1])]  = 1.

    # update lagged variables
    JJ[first(eq[:LR]),first(endo[:R′_t1])] = 1.
    JJ[first(eq[:LR]),first(endo[:R_t])]   = -1.

    JJ[first(eq[:LI]),first(endo[:i′_t1])] = 1.
    JJ[first(eq[:LI]),first(endo[:i_t])]   = -1.

    JJ[first(eq[:LY]),first(endo[:y′_t1])] = 1.
    JJ[first(eq[:LY]),first(endo[:y_t])]   = -1.

    JJ[first(eq[:LW]),first(endo[:w′_t1])] = 1.
    JJ[first(eq[:LW]),first(endo[:w_t])]   = -1.

    JJ[first(eq[:LX]),first(endo[:I′_t1])] = 1.
    JJ[first(eq[:LX]),first(endo[:I_t])]   = -1.

    # discount factor shock
    JJ[first(eq[:eq_b]),first(endo[:b′_t])] = 1.
    JJ[first(eq[:eq_b]),first(endo[:b_t])]  = -ρ_B

    # g/y shock
    JJ[first(eq[:eq_g]),first(endo[:g′_t])] = 1.
    JJ[first(eq[:eq_g]),first(endo[:g_t])]  = -ρ_G

    # tfp growth shock
    JJ[first(eq[:eq_z]),first(endo[:z′_t])] = 1.
    JJ[first(eq[:eq_z]),first(endo[:z_t])]  = -ρ_z

    # investment shock
    JJ[first(eq[:eq_μ]),first(endo[:μ′_t])] = 1.
    JJ[first(eq[:eq_μ]),first(endo[:μ_t])]  = -ρ_μ

    # wage mkup shock
    JJ[first(eq[:eq_λ_w]),first(endo[:λ_w′_t])] = 1.
    JJ[first(eq[:eq_λ_w]),first(endo[:λ_w_t])]  = -ρ_lamw

    # price mkup shock
    JJ[first(eq[:eq_λ_f]),first(endo[:λ_f′_t])] = 1.
    JJ[first(eq[:eq_λ_f]),first(endo[:λ_f_t])]  = -ρ_lamf

    # monetary policy shock
    JJ[first(eq[:eq_rm]),first(endo[:rm′_t])] = 1.
    JJ[first(eq[:eq_rm]),first(endo[:rm_t])]  = -ρ_mon

    #consumption
    #=JJ[first(eq[:eq_consumption]), endo[:l_t]] = (μ .*unc.*xswts.*c)'
    JJ[first(eq[:eq_consumption]), endo[:kf_t]] = -(xswts.*c)' # note, now we linearize
    JJ[first(eq[:eq_consumption]),first(endo[:R_t])]   = -(xswts.*c)'*dF2_dRZ
    JJ[first(eq[:eq_consumption]),first(endo[:z_t])]   = (xswts.*c)'*dF2_dRZ
    JJ[first(eq[:eq_consumption]),first(endo[:w_t])]   = -(xswts.*c)'*dF2_dWH
    JJ[first(eq[:eq_consumption]),first(endo[:L_t])]   = -(xswts.*c)'*dF2_dWH
    JJ[first(eq[:eq_consumption]),first(endo[:t_t])]   = -(xswts.*c)'*dF2_dTT=#

    if !m.testing && get_setting(m, :normalize_distr_variables)
        JJ = normalize(m, JJ)
    end
    return JJ, dF2_dRZ, dF2_dWH, dF2_dTT
end

function euler_equation_hetdsge(nx::Int, ns::Int,
                                qp::Function, qfunction::Function,
                                xgrid::Vector{Float64}, sgrid::Vector{Float64},
                                fgrid::Matrix{Float64},
                                unc::BitArray,
                                xswts::Vector{Float64},
                                R::Float64, γ::Float64, β::Float64,
                                η::Float64, ell::Vector{Float64}, T::Float64,
                                ω::Float64, H::Float64)
    nxns = nx*ns
    ee  = zeros(nxns,nxns) # ee[i,j] takes you from i to j
    ξ   = zeros(nxns,nxns)
    Ξ   = zeros(nxns,nxns)
    dF1_dELL = zeros(nxns,nxns)
    dF1_dRZ  = zeros(nxns)
    dF1_dELLP = zeros(nxns,nxns)
    dF1_dWHP = zeros(nxns)
    dF1_dTTP = zeros(nxns)
    for iss=1:ns
        for ia=1:nx
            i  = nx*(iss-1)+ia
            sumELL = 0.
            sumRZ  = 0.
            sumWH  = 0.
            sumTT  = 0.
            for isp=1:ns
                for iap=1:nx
                    ip = nx*(isp-1)+iap
                    ee[i,ip] = (xgrid[iap] - R*(exp(-γ))*max(xgrid[ia]-1/ell[i], -η) - T)/(ω*H*sgrid[isp])
                    ξ[i,ip] = ((β*R*xswts[i]*exp(-γ))/(ω*H*sgrid[isp])^2)*max(ell[ip],1/(xgrid[iap]+η))*qp(ee[i,ip])*fgrid[iss,isp]
                    Ξ[i,ip] = ((β*R*xswts[i]*exp(-γ))/(ω*H*sgrid[isp]))*max(ell[ip],1/(xgrid[iap]+η))*qfunction(ee[i,ip])*fgrid[iss,isp]
                    sumELL += ξ[i,ip]*R*(exp(-γ))*unc[i]/ell[i]
                    sumRZ  += ξ[i,ip]*R*(exp(-γ))*max(xgrid[ia] - 1/ell[i],-η)
                    dF1_dELLP[i,ip] = Ξ[i,ip]*unc[ip]
                    sumWH  += Ξ[i,ip] + ξ[i,ip]*ee[i,ip]*(ω*H*sgrid[isp])
                    sumTT  += ξ[i,ip]*T
                end
            end
            dF1_dELL[i,i] = -ell[i] - sumELL
            dF1_dRZ[i] = ell[i] - sumRZ
            dF1_dWHP[i] = -sumWH
            dF1_dTTP[i] = -sumTT
        end
    end
    return dF1_dELL, dF1_dRZ, dF1_dELLP, dF1_dWHP, dF1_dTTP, ee
end

function kolmogorov_fwd_hetdsge(nx::Int, ns::Int,
                                qfunction::Function, qp::Function,
                                xgrid::Vector{Float64}, sgrid::Vector{Float64},
                                fgrid::Matrix{Float64}, unc::BitArray,
                                xswts::Vector{Float64},
                                R::Float64, γ::Float64,
                                ell::Vector{Float64}, μ::Vector{Float64},
                                η::Float64, T::Float64, ω::Float64, H::Float64, ee::Matrix{Float64})
    nxns = nx*ns
    bigΨ    = zeros(nxns,nxns)
    smallψ = zeros(nxns,nxns)
    dF2_dRZ = zeros(nxns)
    dF2_dM = zeros(nxns,nxns)
    dF2_dELL = zeros(nxns,nxns)
    dF2_dWH = zeros(nxns)
    dF2_dTT = zeros(nxns)

    for isp=1:ns
        for iap=1:nx
            ip  = nx*(isp-1)+iap
            sumWH = 0.
            sumRZ = 0.
            sumTT = 0.
            for iss=1:ns
                for ia=1:nx
                    i = nx*(iss-1)+ia
                    bigΨ[ip,i] = xswts[i]*μ[i]*qfunction(ee[i,ip])*fgrid[iss,isp]/(ω*H*sgrid[isp])
                    smallψ[ip,i] = xswts[i]*μ[i]*qp(ee[i,ip])*fgrid[iss,isp]/((ω*H*sgrid[isp])^2)
                    sumWH += bigΨ[ip,i] + smallψ[ip,i]*ee[i,ip]*(ω*H*sgrid[isp])
                    sumRZ += smallψ[ip,i]*(R*exp(-γ))*max(xgrid[ia] - 1/ell[i],-η)
                    dF2_dELL[ip,i] = smallψ[ip,i]*(R*exp(-γ))*(unc[i]/ell[i])
                    sumTT += smallψ[ip,i]*T
                    dF2_dM[ip,i] = xswts[i]*qfunction(ee[i,ip])*fgrid[iss,isp]/(ω*H*sgrid[isp]) # note, now we linearize
                end
            end
            dF2_dWH[ip] = -sumWH
            dF2_dRZ[ip] = -sumRZ
            dF2_dTT[ip] = -sumTT
        end
    end
    return dF2_dWH, dF2_dRZ, dF2_dTT, dF2_dELL, bigΨ, dF2_dM
end

function normalize(m::HetDSGE, JJ::Matrix{Float64})

    Qx, _, Qleft, Qright = compose_normalization_matrices(m)

    m <= Setting(:n_predetermined_variables, size(Qx, 1))

	Jac1 = Qleft*JJ*Qright

    return Jac1
end

function compose_normalization_matrices(m::HetDSGE)
    nx::Int = get_setting(m, :nx)
    ns::Int = get_setting(m, :ns)
    nscalars::Int = get_setting(m, :nscalars)
    nyscalars::Int = get_setting(m, :nyscalars)
    nxscalars::Int = get_setting(m, :nxscalars)

    # Create PPP matrix
    P1 = kron(Matrix{Float64}(I, ns,ns),ones(nx,1))
    Ptemp = Matrix{Float64}(I, nx, nx)
    Ptemp = Ptemp[:, 2:end]
    P2 = kron(Matrix{Float64}(I, ns, ns), Ptemp)
    P  = hcat(P1, P2)

    Q,R = qr(P)
    Q = Array(Q)
    S         = Q[:, ns+1:end]'

    nxns = nx*ns

    Qleft     = cat(Matrix{Float64}(I, nxns, nxns),S,Matrix{Float64}(I, nscalars, nscalars), dims = [1 2])
    Qx        = cat(S,Matrix{Float64}(I, nxscalars, nxscalars), dims = [1 2])
    Qy        = cat(Matrix{Float64}(I, nxns, nxns),Matrix{Float64}(I, nyscalars, nyscalars), dims = [1 2])
    Qright    = cat(Qx',Qy',Qx',Qy', dims = [1,2])

    return Qx, Qy, Qleft, Qright
end

function truncate_distribution!(m::HetDSGE)
    mindens = get_setting(m, :mindens)
    trunc_distr = get_setting(m, :trunc_distr)
    rescale_weights = get_setting(m, :rescale_weights)

    nx = get_setting(m, :nx)
    μ = m[:μstar].value
    ell = m[:lstar].value
    c = m[:cstar].value
    xgrid = m.grids[:xgrid].points
    xlo = get_setting(m, :xlo)
    swts::Vector{Float64}  = m.grids[:sgrid].weights

    if trunc_distr
        oldnx = nx
        nx = maximum(findall(μ[1:nx]+μ[nx+1:2*nx] .> mindens)) # used to be 1e-8
        m[:μstar] = μ[[1:nx;oldnx+1:oldnx+nx]]
        m[:lstar] = ell[[1:nx;oldnx+1:oldnx+nx]]
        m[:cstar] = c[[1:nx;oldnx+1:oldnx+nx]]
        if rescale_weights
            xhi = xgrid[nx]
            xscale = xhi-xlo
        end
        m <= Setting(:nx, nx)
        m <= Setting(:xhi, xhi)
        m <= Setting(:xscale, xscale)
        m.grids[:xgrid] = Grid(uniform_quadrature(xscale), xlo, xhi, nx, scale = xscale)
        m.grids[:weights_total] = kron(swts, m.grids[:xgrid].weights)
        nxns = nx*get_setting(m, :ns)
        m <= Setting(:n, nxns)
        m <= Setting(:nx, nx)

        setup_indices!(m)

        init_states_and_jumps!(m, get_setting(m, :states), get_setting(m, :jumps))

        normalize_model_state_indices!(m)

        endogenous_states_augmented = [:i_t1, :c_t, :c_t1]
        for (i,k) in enumerate(endogenous_states_augmented); m.endogenous_states_augmented[k] = i + first(collect(values(m.endogenous_states))[end]) end


    end
end
