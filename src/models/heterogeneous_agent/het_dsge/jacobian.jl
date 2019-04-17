function jacobian(m::HetDSGE)

    # Load in endogenous state and eq cond indices
    endo = m.endogenous_states #=augment_model_states(m.endogenous_states_unnormalized,
                                n_model_states_unnormalized(m)) =#
    eq   = m.equilibrium_conditions

    # Load in parameters, steady-state parameters, and grids
    r::Float64     = m[:r].value
    α::Float64  = m[:α].value
    H::Float64     = m[:H].value
    δ::Float64     = m[:δ].value
    μ_sp::Float64    = m[:μ_sp].value
    ρ_sp::Float64  = m[:ρ_sp].value
    σ_sp::Float64  = m[:σ_sp].value
    γ::Float64  = m[:γ].value
    g::Float64     = m[:g].value
    η::Float64 = m[:η].value
    ρB::Float64 = m[:ρB].value
    ρG::Float64 = m[:ρG].value
    ρZ::Float64 = m[:ρZ].value
    ρμ::Float64 = m[:ρμ].value
    ρlamw::Float64 = m[:ρlamw].value
    ρlamf::Float64 = m[:ρlamf].value
    ρmon::Float64 = m[:ρmon].value
    spp::Float64 = m[:spp].value
    lamw::Float64 = m[:lamw].value
    ϕh::Float64 = m[:ϕh].value
    Φw::Float64 = m[:Φw].value
    lamf::Float64 = m[:lamf].value
    Φp::Float64 = m[:Φp].value
    ρR::Float64 = m[:ρR].value
    ψπ::Float64 = m[:ψπ].value
    ψy::Float64 = m[:ψy].value

    R = 1 + r

    ell::Vector{Float64}  = m[:lstar].value
    c::Vector{Float64}    = m[:cstar].value
    μ::Vector{Float64}    = m[:μstar].value
    β::Float64            = m[:βstar].value

    T::Float64 = m[:Tstar].value
    ω::Float64 = m[:ωstar].value
    xstar::Float64 = m[:xstar].value
    ystar::Float64 = m[:ystar].value
    Rk::Float64 = m[:Rkstar].value
    kstar::Float64 = m[:kstar].value


    xgrid::Vector{Float64} = m.grids[:xgrid].points
    xwts::Vector{Float64}  = m.grids[:xgrid].weights
    sgrid::Vector{Float64} = m.grids[:sgrid].points
    swts::Vector{Float64}  = m.grids[:sgrid].weights
    fgrid::Matrix{Float64} = m.grids[:fgrid]
    xgrid_total::Vector{Float64} = m.grids[:xgrid_total]
    sgrid_total::Vector{Float64} = m.grids[:sgrid_total]
    weights_total::Vector{Float64} = m.grids[:weights_total]

    xswts = kron(swts,xwts)

    zlo::Float64 = get_setting(m, :zlo)
    zhi::Float64 = get_setting(m, :zhi)

    nx::Int = get_setting(m, :nx)
    ns::Int = get_setting(m, :ns)
    nxns = nx*ns

    zgrid = range(zlo,stop = zhi,length = nx)
    zwts = (zhi-zlo)/nx
    sumz=0.
    for i=1:nx
        sumz += mollifier_hetdsge(zgrid[i],zhi,zlo)*zwts
    end

    qp(z) = dmollifier_hetdsge(z, zhi, zlo)
    qfunction(x) = mollifier_hetdsge(x, zhi, zlo)/sumz

    unc = 1 ./ ell .<= repeat(xgrid,ns) .+ η

    # Euler Equation
    dF1_dELL, dF1_dRZ, dF1_dELLP, dF1_dWHP, dF1_dTTP, ee = euler_equation_hetdsge(nx, ns, qp, qfunction, xgrid, sgrid, sgrid_total, fgrid, unc,
                           xwts, swts, xswts,
                           R, γ, β,
                           c, η, ell, T, ω, H)

    # KF Equation
    dF2_dWH, dF2_dRZ, dF2_dTT,dF2_dELL, bigΨ = kolmogorov_fwd_hetdsge(nx, ns, qfunction, qp, xgrid, sgrid, fgrid, unc, xwts, swts, xswts, R, γ,
                       ell, μ, η, c, T, ω, H, ee)

    # Market clearing, lambda function
    c = min.(1 ./ ell,repeat(xgrid,ns).+η)
    lam = (xswts.*μ)'*(1./c) # average marginal utility which the union uses to set wages
    ϕ = lam*ω/(H^ϕh) # now that we know lam in steady state, choose disutility to target hours H


    nvars = get_setting(m, :nvars)
    # Make the Jacobian
    JJ = zeros(nvars, 2*nvars)

    # Euler equation
    JJ[eq[:eq_euler],endo[:l′_t]] = dF1_dELLP
    JJ[eq[:eq_euler],endo[:ZP]]   = -dF1_dRZ
    JJ[eq[:eq_euler],endo[:w′_t]]   = dF1_dWHP
    JJ[eq[:eq_euler],endo[:L′_t]]  = dF1_dWHP
    JJ[eq[:eq_euler],endo[:t′_t]]  = dF1_dTTP
    JJ[eq[:eq_euler],endo[:B]]    = ell
    JJ[eq[:eq_euler],endo[:ELL]]  = dF1_dELL
    JJ[eq[:eq_euler],endo[:R_t]]   = dF1_dRZ

    # KF eqn
    JJ[eq[:eq_kolmogorov_fwd],endo[:μ_t1]]   = bigΨ
    JJ[eq[:eq_kolmogorov_fwd],endo[:R_t1]]  = dF2_dRZ
    JJ[eq[:eq_kolmogorov_fwd],endo[:l_t1]] = dF2_dELL
    JJ[eq[:eq_kolmogorov_fwd],endo[:Z]]    = -dF2_dRZ
    JJ[eq[:eq_kolmogorov_fwd],endo[:μ_t]]    = -Matrix(Diagonal(μ))
    JJ[eq[:eq_kolmogorov_fwd],endo[:w_t]]    = dF2_dWH
    JJ[eq[:eq_kolmogorov_fwd],endo[:L_t]]   = dF2_dWH
    JJ[eq[:eq_kolmogorov_fwd],endo[:t_t]]   = dF2_dTT

    # update lagged ELL
    JJ[eq[:eq_lag_ell],endo[:l′_t1]] = eye(nxns)
    JJ[eq[:eq_lag_ell],endo[:ELL]]   = -eye(nxns)

    # update lagged M
    JJ[eq[:eq_lag_wealth],endo[:μ′_t1]] = eye(nxns)
    JJ[eq[:eq_lag_wealth],endo[:μ_t]]   = -eye(nxns)

    # mkt clearing
    JJ[first(eq[:eq_market_clearing]),first(endo[:y_t])]   = ystar/g
    JJ[first(eq[:eq_market_clearing]),first(endo[:G])]   = -ystar/g
    JJ[first(eq[:eq_market_clearing]),first(endo[:I_t])]   = -xstar
    JJ[eq[:eq_market_clearing],endo[:ELL]] = (μ.*unc.*xswts.*c)'
    JJ[eq[:eq_market_clearing],endo[:μ_t]]   = -(μ.*xswts.*c)'

    # lambda = average marginal utility
    JJ[first(eq[:eq_lambda]),first(endo[:mu_t])] = lam
    JJ[first(eq[:eq_lambda]),endo[:μ_t]]   = -(xswts.*μ./c)'
    JJ[first(eq[:eq_lambda]),endo[:ELL]] = -(xswts.*unc.*μ./c)'

    # transfer
    JJ[first(eq[:eq_transfers]),first(endo[:t_t])]  = T
    JJ[first(eq[:eq_transfers]),first(endo[:capreturn_t])]  = -Rk*kstar
    JJ[first(eq[:eq_transfers]),first(endo[:k_t])]  = -Rk*kstar
    JJ[first(eq[:eq_transfers]),first(endo[:Z])]  = Rk*kstar
    JJ[first(eq[:eq_transfers]),first(endo[:I_t])]   = xstar
    JJ[first(eq[:eq_transfers]),first(endo[:mc_t])]  = ystar
    JJ[first(eq[:eq_transfers]),first(endo[:y_t])]   = (1-1/g)
    JJ[first(eq[:eq_transfers]),first(endo[:G])]   = (ystar/g)

    # investment
    JJ[first(eq[:eq_investment]),first(endo[:Q_t])]  = 1.
    JJ[first(eq[:eq_investment]),first(endo[:MU])] = 1.
    JJ[first(eq[:eq_investment]),first(endo[:I′_t])] = spp*(exp(3*γ))/R
    JJ[first(eq[:eq_investment]),first(endo[:ZP])] = spp*(exp(3*γ))/R
    JJ[first(eq[:eq_investment]),first(endo[:I_t])]  = -spp*(exp(3*γ))/R - spp*exp(2*γ)
    JJ[first(eq[:eq_investment]),first(endo[:I_t1])] = spp*exp(2*γ)
    JJ[first(eq[:eq_investment]),first(endo[:Z])]  = -spp*exp(2*γ)

    # tobin's q
    JJ[first(eq[:eq_tobin_q]),first(endo[:R_t])]  = R
    JJ[first(eq[:eq_tobin_q]),first(endo[:Q_t])]   = R
    JJ[first(eq[:eq_tobin_q]),first(endo[:capreturn′_t])] = -Rk
    JJ[first(eq[:eq_tobin_q]),first(endo[:Q′_t])]  = -(1-δ)

    # capital accumulation
    JJ[first(eq[:eq_capital_accumulation]),first(endo[:k′_t])] = 1.
    JJ[first(eq[:eq_capital_accumulation]),first(endo[:k_t])]  = -(1-δ)
    JJ[first(eq[:eq_capital_accumulation]),first(endo[:Z])]   = (1-δ)
    JJ[first(eq[:eq_capital_accumulation]),first(endo[:MU])]  = -xstar/kstar
    JJ[first(eq[:eq_capital_accumulation]),first(endo[:I_t])]   = -xstar/kstar

    # wage phillips curve
    JJ[first(eq[:eq_wage_phillips]),first(endo[:wageinflation_t])]  = -1.
    JJ[first(eq[:eq_wage_phillips]),first(endo[:LAMW])] = (ϕ*H^ϕh)/Φw
    JJ[first(eq[:eq_wage_phillips]),first(endo[:L_t])]   = (ϕ*(H^ϕh)*(1+lamw)/lamw*Φw)*ϕh
    JJ[first(eq[:eq_wage_phillips]),first(endo[:mu_t])]  = -(ϕ*(H^ϕh)*(1+lamw)/lamw*Φw)
    JJ[first(eq[:eq_wage_phillips]),first(endo[:w_t])]    = -(ϕ*(H^ϕh)*(1+lamw)/lamw*Φw)
    JJ[first(eq[:eq_wage_phillips]),first(endo[:wageinflation′_t])]  = β

    # price phillips curve
    JJ[first(eq[:eq_price_phillips]),first(endo[:π_t])]   = -1.
    JJ[first(eq[:eq_price_phillips]),first(endo[:mc_t])]   =(1+lamf)/(lamf*Φp)
    JJ[first(eq[:eq_price_phillips]),first(endo[:LAMF])] = 1/Φp
    JJ[first(eq[:eq_price_phillips]),first(endo[:π′_t])]  = 1/R

    # marginal cost
    JJ[first(eq[:eq_marginal_cost]),first(endo[:mc_t])] = 1.
    JJ[first(eq[:eq_marginal_cost]),first(endo[:w_t])]  = -(1-α)
    JJ[first(eq[:eq_marginal_cost]),first(endo[:capreturn_t])] = -α

    # gdp
    JJ[first(eq[:eq_gdp]),first(endo[:y_t])]  = 1.
    JJ[first(eq[:eq_gdp]),first(endo[:Z])]  = α
    JJ[first(eq[:eq_gdp]),first(endo[:k_t])] = -α
    JJ[first(eq[:eq_gdp]),first(endo[:L_t])] = -(1-α)

    # optimal k/l ratio
    JJ[first(eq[:eq_optimal_kl]),first(endo[:capreturn_t])] = 1.
    JJ[first(eq[:eq_optimal_kl]),first(endo[:w_t])]  = -1.
    JJ[first(eq[:eq_optimal_kl]),first(endo[:L_t])] = -1.
    JJ[first(eq[:eq_optimal_kl]),first(endo[:k_t])] = 1.
    JJ[first(eq[:eq_optimal_kl]),first(endo[:Z])]  = -1.

    # taylor rule
    JJ[first(eq[:eq_taylor]),first(endo[:R_t1])]   = -1.
    JJ[first(eq[:eq_taylor]),first(endo[:i_t1])]  = ρR
    JJ[first(eq[:eq_taylor]),first(endo[:π_t])]   = (1-ρR)*ψπ
    JJ[first(eq[:eq_taylor]),first(endo[:π_t1])]  = (1-ρR)*ψπ
    JJ[first(eq[:eq_taylor]),first(endo[:π_t2])] = (1-ρR)*ψπ
    JJ[first(eq[:eq_taylor]),first(endo[:π_t3])] = (1-ρR)*ψπ
    JJ[first(eq[:eq_taylor]),first(endo[:y_t])]    = (1-ρR)*ψy
    JJ[first(eq[:eq_taylor]),first(endo[:y_t4])]  = -(1-ρR)*ψy
    JJ[first(eq[:eq_taylor]),first(endo[:Z])]    = (1-ρR)*ψy
    JJ[first(eq[:eq_taylor]),first(endo[:z_t1])]   = (1-ρR)*ψy
    JJ[first(eq[:eq_taylor]),first(endo[:z_t2])]  = (1-ρR)*ψy
    JJ[first(eq[:eq_taylor]),first(endo[:z_t3])]  = (1-ρR)*ψy
    JJ[first(eq[:eq_taylor]),first(endo[:MON])]  = 1.

    # fisher eqn
    JJ[first(eq[:eq_fisher]),first(endo[:R_t])]  = 1.
    JJ[first(eq[:eq_fisher]),first(endo[:π′_t])] = 1.
    JJ[first(eq[:eq_fisher]),first(endo[:I_t])]  = -1.

    # wage inflation
    JJ[first(eq[:eq_nominal_wage_inflation]),first(endo[:wageinflation_t])] = 1.
    JJ[first(eq[:eq_nominal_wage_inflation]),first(endo[:π_t])]  = -1.
    JJ[first(eq[:eq_nominal_wage_inflation]),first(endo[:Z])]   = -1.
    JJ[first(eq[:eq_nominal_wage_inflation]),first(endo[:w_t])]   = -1.
    JJ[first(eq[:eq_nominal_wage_inflation]),first(endo[:w_t1])]  = 1.

    # update lagged variables
    JJ[first(eq[:LR]),first(endo[:R′_t1])] = 1.
    JJ[first(eq[:LR]),first(endo[:R_t])]   = -1.

    JJ[first(eq[:LI]),first(endo[:i′_t1])] = 1.
    JJ[first(eq[:LI]),first(endo[:i_t])]   = -1.

    JJ[first(eq[:LPI]),first(endo[:π′_t1])] = 1.
    JJ[first(eq[:LPI]),first(endo[:π_t])]   = -1.

    JJ[first(eq[:L2PI]),first(endo[:π′_t2])] = 1.
    JJ[first(eq[:L2PI]),first(endo[:π_t1])]   = -1.

    JJ[first(eq[:L3PI]),first(endo[:π′_t3])] = 1.
    JJ[first(eq[:L3PI]),first(endo[:π_t2])]  = -1.

    JJ[first(eq[:LY]),first(endo[:y′_t1])] = 1.
    JJ[first(eq[:LY]),first(endo[:y_t])]   = -1.

    JJ[first(eq[:L2Y]),first(endo[:y′_t2])] = 1.
    JJ[first(eq[:L2Y]),first(endo[:y_t1])]   = -1.

    JJ[first(eq[:L3Y]),first(endo[:y′_t3])] = 1.
    JJ[first(eq[:L3Y]),first(endo[:y_t2])]  = -1.

    JJ[first(eq[:L4Y]),first(endo[:y′_t4])] = 1.
    JJ[first(eq[:L4Y]),first(endo[:y_t3])]  = -1.

    JJ[first(eq[:LZ]),first(endo[:z′_t1])] = 1.
    JJ[first(eq[:LZ]),first(endo[:Z])]   = -1.

    JJ[first(eq[:L2Z]),first(endo[:z′_t2])] = 1.
    JJ[first(eq[:L2Z]),first(endo[:z_t1])]   = -1.

    JJ[first(eq[:L3Z]),first(endo[:z′_t3])] = 1.
    JJ[first(eq[:L3Z]),first(endo[:z_t2])]  = -1.

    JJ[first(eq[:LW]),first(endo[:w′_t1])] = 1.
    JJ[first(eq[:LW]),first(endo[:w_t])]   = -1.

    JJ[first(eq[:LX]),first(endo[:I′_t1])] = 1.
    JJ[first(eq[:LX]),first(endo[:I_t])]   = -1.

    # discount factor shock
    JJ[first(eq[:F33]),first(endo[:BP])] = 1.
    JJ[first(eq[:F33]),first(endo[:B])]  = -ρB

    # g/y shock
    JJ[first(eq[:F34]),first(endo[:GP])] = 1.
    JJ[first(eq[:F34]),first(endo[:G])]  = -ρG

    # tfp growth shock
    JJ[first(eq[:F35]),first(endo[:ZP])] = 1.
    JJ[first(eq[:F35]),first(endo[:Z])]  = -ρZ

    # investment shock
    JJ[first(eq[:F36]),first(endo[:MUP])] = 1.
    JJ[first(eq[:F36]),first(endo[:MU])]  = -ρμ

    # wage mkup shock
    JJ[first(eq[:F37]),first(endo[:LAMWP])] = 1.
    JJ[first(eq[:F37]),first(endo[:LAMW])]  = -ρlamw

    # price mkup shock
    JJ[first(eq[:F38]),first(endo[:LAMFP])] = 1.
    JJ[first(eq[:F38]),first(endo[:LAMF])]  = -ρlamf

    # monetary policy shock
    JJ[first(eq[:F39]),first(endo[:MONP])] = 1.
    JJ[first(eq[:F39]),first(endo[:MON])]  = -ρmon

    if !m.testing && get_setting(m, :normalize_distr_variables)
        JJ = normalize(m, JJ)
    end

    return JJ
end

function euler_equation_hetdsge(nx::Int, ns::Int,
                                qp::Function, qfunction::Function,
                                xgrid::Vector{Float64}, sgrid::Vector{Float64},
                                sgrid_total::Vector{Float64}, fgrid::Matrix{Float64},
                                unc::BitArray,
                                xwts::Vector{Float64}, swts::Vector{Float64}, xswts::Vector{Float64},
                                R::Float64, γ::Float64, β::Float64,
                                c::Vector{Float64}, η::Float64, ell::Vector{Float64}, T::Float64,
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
                    sumWH  += Ξ[i,ip] + ξ[i,ip]*ee[i,ip]/(ω*H*sgrid[isp])
                    sumTT  += ξ[i,ip]*T
                end
            end
            dF1_dELL[i,i] = -ell[i] - sumELL
            dF1_dRZ[i] = ell[i] - sumRZ
            dF1_dWHP[i] = -sumWH
            dF1_dTTP[i] = sumTT
        end
    end
    return dF1_dELL, dF1_dRZ, dF1_dELLP, dF1_dWHP, dF1_dTTP, ee
end

function kolmogorov_fwd_hetdsge(nx::Int, ns::Int,
                                qfunction::Function, qp::Function,
                                xgrid::Vector{Float64}, sgrid::Vector{Float64}, fgrid::Matrix{Float64}, unc::BitArray,
                                xwts::Vector{Float64}, swts::Vector{Float64}, xswts::Vector{Float64},
                                R::Float64, γ::Float64,
                                ell::Vector{Float64}, μ::Vector{Float64}, η::Float64, c::Vector{Float64}, T::Float64, ω::Float64, H::Float64, ee::Matrix{Float64})
    nxns = nx*ns
    bigΨ    = zeros(nxns,nxns) # this is also dF2_dM
    smallψ = zeros(nxns,nxns)
    dF2_dRZ = zeros(nxns)
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
                    sumWH += bigΨ[ip,i] + smallψ[ip,i]*ee[i,ip]/(ω*H*sgrid[isp])
                    sumRZ += smallψ[ip,i]*(R*exp(-γ))*max(xgrid[ia] - 1/ell[i],-η)
                    dF2_dELL[ip,i] = smallψ[ip,i]*(R*exp(-γ))*(unc[i]/ell[i])
                    sumTT += smallψ[ip,i]*T
                end
            end
            dF2_dWH[ip] = -sumWH
            dF2_dRZ[ip] = -sumRZ
            dF2_dTT[ip] = -sumTT
        end
    end
    return dF2_dWH, dF2_dRZ, dF2_dTT, dF2_dELL, bigΨ
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
    P1 = kron(eye(ns),ones(nx,1))
    Ptemp = eye(nx)
    Ptemp = Ptemp[:, 2:end]
    P2 = kron(eye(ns), Ptemp)
    P  = hcat(P1, P2)

    Q,R = qr(P)
    Q = Array(Q)
    S         = Q[:, ns+1:end]'

    nxns = nx*ns

    Qleft     = cat(eye(nxns),S,eye(nxns),S,eye(nscalars), dims = [1 2])
    Qx        = cat(eye(nxns),S,eye(nxscalars), dims = [1 2])
    Qy        = cat(eye(nxns),S,eye(nyscalars), dims = [1 2])
    Qright    = cat(Qx',Qy',Qx',Qy', dims = [1,2])

    return Qx, Qy, Qleft, Qright
end
