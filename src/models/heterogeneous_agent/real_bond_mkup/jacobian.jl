function jacobian(m::RealBondMkup)

    # Load in endogenous state and eq cond indices
    endo = augment_model_states(m.endogenous_states_unnormalized,
                                n_model_states_unnormalized(m))
    eq   = m.equilibrium_conditions

    # Load in parameters, steady-state parameters, and grids
    γ::Float64       = m[:γ].value
    abar::Float64    = m[:abar].value
    R::Float64       = m[:R].value
    ν::Float64       = m[:ν].value
    ρ_z::Float64     = m[:ρ_z].value
    ρ_mon::Float64   = m[:ρ_mon].value
    ρ_mkp::Float64   = m[:ρ_mkp].value
    ρ_tay::Float64   = m[:ρ_tay].value
    κ::Float64       = m[:κ].value
    phipi::Float64   = m[:phipi].value
    aborrow::Float64 = abar / R

    ell::Vector{Float64}  = m[:lstar].value
    c::Vector{Float64}    = m[:cstar].value
    μ::Vector{Float64}    = m[:μstar].value
    η::Vector{Float64}    = m[:ηstar].value
    χss::Vector{Float64}  = m[:χstar].value
    β::Float64            = m[:βstar].value

    xgrid::Vector{Float64} = m.grids[:xgrid].points
    xwts::Vector{Float64}  = m.grids[:xgrid].weights
    sgrid::Vector{Float64} = m.grids[:sgrid].points
    swts::Vector{Float64}  = m.grids[:sgrid].weights
    ggrid::Vector{Float64} = m.grids[:ggrid]
    xgrid_total::Vector{Float64} = m.grids[:xgrid_total]
    sgrid_total::Vector{Float64} = m.grids[:sgrid_total]
    weights_total::Vector{Float64} = m.grids[:weights_total]

    elo::Float64 = get_setting(m, :elo)
    ehi::Float64 = get_setting(m, :ehi)

    nx::Int = get_setting(m, :nx)
    ns::Int = get_setting(m, :ns)

    # Mollifiers
    qfunction(x) = mollifier_realbond(x,ehi,elo)
    qp(x) = dmollifier_realbond(x,ehi,elo)

    # Jacobian stuff
    chipW, chipR, chipX = construct_chip_realbondmkup(xgrid_total, γ, ν, aborrow, abar, R, χss)

    # Create auxiliary variables
    # c = min.(ell.^(-1.0/γ),χss) # Consumption Decision Rule
    # η = (sgrid_total.^(1.0/ν)).*(c.^(-γ/ν))

    # EE
    dF1_dELLP, dF1_dRRP, dF1_dWWP, dF1_dTTP, dF1_dELL, dF1_dWW, dF1_dRR, dF1_dTT, unc, a =
    euler_equation_realbondmkup(nx, ns, qp, qfunction, xgrid, sgrid, ggrid, sgrid_total, xwts, swts,
                            chipR, chipW, chipX, R, γ, β, ν,
                            c, η, ell, χss)

    # KF
    dF2_dMU, dF2_dELL, dF2_dRR, dF2_dWW, dF2_dTT =
    kolmogorov_fwd_realbond(nx, ns, qfunction, qp, xgrid, sgrid, ggrid, xwts, swts,
                            unc, a, R, ν, γ, ell, μ, η, c)

    # MKT clearing
    dF3_dMU, dF3_dELL, dF3_dWW, dF3_dRR, dF3_dTT, GDP =
    market_clearing_realbondmkup(nx, ns, xgrid, sgrid, xwts, swts,
                             η, c, μ, ell, γ, ν, chipW, chipR, chipX, unc)

    # Make the Jacobian
    JJ = zeros(2*nx*ns+9,4*nx*ns+18)

    # Euler equation
    JJ[eq[:eq_euler], endo[:l′_t]]         = dF1_dELLP
    JJ[eq[:eq_euler], first(endo[:R′_t])]  = dF1_dRRP
    JJ[eq[:eq_euler], first(endo[:w′_t])]  = dF1_dWWP
    JJ[eq[:eq_euler], first(endo[:t′_t])]  = dF1_dTTP
    JJ[eq[:eq_euler], endo[:l_t]]          = dF1_dELL
    JJ[eq[:eq_euler], first(endo[:R_t])]   = dF1_dRR
    JJ[eq[:eq_euler], first(endo[:w_t])]   = dF1_dWW
    JJ[eq[:eq_euler], first(endo[:t_t])]   = dF1_dTT

    # KF
    JJ[eq[:eq_kolmogorov_fwd], endo[:μ′_t]]        = -eye(nx*ns)
    JJ[eq[:eq_kolmogorov_fwd], endo[:μ_t]]         = dF2_dMU
    JJ[eq[:eq_kolmogorov_fwd], endo[:l_t]]         = dF2_dELL
    JJ[eq[:eq_kolmogorov_fwd], first(endo[:w_t])]  = dF2_dWW
    JJ[eq[:eq_kolmogorov_fwd], first(endo[:R_t])]  = dF2_dRR
    JJ[eq[:eq_kolmogorov_fwd], first(endo[:t_t])]  = dF2_dTT

    # mkt ckr
    JJ[first(eq[:eq_market_clearing]), endo[:μ_t]]         = dF3_dMU
    JJ[first(eq[:eq_market_clearing]), first(endo[:z_t])]  = GDP
    JJ[first(eq[:eq_market_clearing]), endo[:l_t]]         = dF3_dELL
    JJ[first(eq[:eq_market_clearing]), first(endo[:w_t])]  = dF3_dWW
    JJ[first(eq[:eq_market_clearing]), first(endo[:R_t])]  = dF3_dRR
    JJ[first(eq[:eq_market_clearing]), first(endo[:t_t])]  = dF3_dTT

    # TFP
    JJ[first(eq[:eq_TFP]), first(endo[:z′_t])]  = 1.0
    JJ[first(eq[:eq_TFP]), first(endo[:z_t])]   = -ρ_z

    # Phillips
    JJ[first(eq[:eq_phillips]), first(endo[:π′_t])]  = -1.0/R
    JJ[first(eq[:eq_phillips]), first(endo[:z_t])]   = κ
    JJ[first(eq[:eq_phillips]), first(endo[:w_t])]   = -κ
    JJ[first(eq[:eq_phillips]), first(endo[:π_t])]   = 1.0
    JJ[first(eq[:eq_phillips]), first(endo[:mkp_t])] = -1.0

    # taylor
    JJ[first(eq[:eq_taylor]), first(endo[:i_t])]   = 1.0
    JJ[first(eq[:eq_taylor]), first(endo[:π_t])]   = -phipi*R*(1-ρ_tay)
    JJ[first(eq[:eq_taylor]), first(endo[:mon_t])] = -1.0
    JJ[first(eq[:eq_taylor]), first(endo[:l_i_t])] = -ρ_tay

    # fisher
    JJ[first(eq[:eq_fisher]), first(endo[:i_t])]  = 1.0
    JJ[first(eq[:eq_fisher]), first(endo[:π′_t])] = -R
    JJ[first(eq[:eq_fisher]), first(endo[:R_t])]  = -1.0

    # transfers
    JJ[first(eq[:eq_transfers]), first(endo[:t_t])]   = 1.0
    JJ[first(eq[:eq_transfers]), first(endo[:z_t])]   = -GDP
    JJ[first(eq[:eq_transfers]), first(endo[:w_t])]   = GDP

    # MP
    JJ[first(eq[:eq_monetary_policy]), first(endo[:mon′_t])] = 1.0
    JJ[first(eq[:eq_monetary_policy]), first(endo[:mon_t])]  = -ρ_mon

    # Markup
    JJ[first(eq[:eq_markup]), first(endo[:mkp′_t])] = 1.0
    JJ[first(eq[:eq_markup]), first(endo[:mkp_t])]  = -ρ_mkp

    # Lagged Monetary Policy
    JJ[first(eq[:eq_lagged_nominal_rate]), first(endo[:l_i′_t])] = 1.0
    JJ[first(eq[:eq_lagged_nominal_rate]), first(endo[:i_t])]    = -1.0

    if !m.testing && get_setting(m, :normalize_distr_variables)
        JJ = normalize(m, JJ)
    end

    return JJ
end

function construct_chip_realbondmkup(xgrid_total::Vector{Float64},
                                 γ::Float64, ν::Float64, aborrow::Float64,
                                 abar::Float64, R::Float64,
                                 χss::Vector{Float64})
    chipW = (1.0+ν)./(ν./(aborrow .+ χss - xgrid_total) + γ./χss)
    chipR = (1.0/(R^2))*(ν*abar./(aborrow .+ χss - xgrid_total))./(ν./(aborrow .+ χss - xgrid_total) + γ./χss)
    chipX = (ν*χss)./(ν*χss + γ*(aborrow .+ χss - xgrid_total))

    return chipW, chipR, chipX
end

function euler_equation_realbondmkup(nx::Int, ns::Int,
                                 qp::Function, qfunction::Function,
                                 xgrid::Vector{Float64}, sgrid::Vector{Float64},
                                 ggrid::Vector{Float64}, sgrid_total::Vector{Float64},
                                 xwts::Vector{Float64}, swts::Vector{Float64},
                                 chipR::Vector{Float64}, chipW::Vector{Float64}, chipX::Vector{Float64},
                                 R::Float64, γ::Float64, β::Float64, ν::Float64,
                                 c::Vector{Float64}, η::Vector{Float64}, ell::Vector{Float64}, χss::Vector{Float64})
    ξ   = zeros(nx*ns,nx*ns)
    Γ   = zeros(nx*ns)
    Ξ   = zeros(nx*ns,nx*ns)
    dF1_dELLP = zeros(nx*ns,nx*ns)
    dF1_dRRP  = zeros(nx*ns)
    dF1_dWWP  = zeros(nx*ns)
    dF1_dTTP  = zeros(nx*ns)
    dF1_dELL  = zeros(nx*ns,nx*ns)
    dF1_dWW   = zeros(nx*ns)
    dF1_dRR   = zeros(nx*ns)
    dF1_dTT   = zeros(nx*ns)
    unc = zeros(nx*ns)
    a   = zeros(nx*ns)

    # EE jacobian terms
    for ix=1:nx
        for is =1:ns
            i = ix + nx*(is-1)
            sumΓ   = 0.0
            sumRRP = 0.0
            sumWWP = 0.0
            sumTTP = 0.0
            sumWW  = 0.0
            sumTT  = 0.0
            sumξ   = 0.0
            sumELL = 0.0
            a[i] = R*(sgrid[is]*η[i]+xgrid[ix] - c[i])
            for ixp=1:nx
                for isp = 1:ns
                    ip = ixp + nx*(isp-1)
                    unc[ip] = ((ell[ip]^(-1/γ))<=χss[ip]) # =1 if unconstrained
                    mf = R*(sgrid[is]*η[i]+xgrid[ix] - c[i]) - xgrid[ixp]
                    ξ[i,ip] = β*R*(c[ip]^(-γ))*qp(mf)*ggrid[isp]*swts[isp]*xwts[ixp]
                    Ξ[i,ip] = β*R*(c[ip]^(-γ-1))*qfunction(mf)*ggrid[isp]*swts[isp]*xwts[ixp]
                    dF1_dELLP[i,ip] = Ξ[i,ip]*(ell[ip]^(-1- 1/γ))*unc[ip]
                    sumRRP += -Ξ[i,ip]*chipR[ip]*unc[ip]
                    sumWWP += -Ξ[i,ip]*chipW[ip]*unc[ip]
                    sumTTP += -Ξ[i,ip]*chipX[ip]*unc[ip]
                    sumELL += ξ[i,ip]*R*unc[i]*((1/ν) + sgrid[is]*η[i]/(ν*c[i]))*ell[i]^(-(1/γ) - 1)
                    sumWW += ξ[i,ip]*R*unc[i]*(1.0+1.0/ν)*sgrid[is]*η[i]
                    sumTT += ξ[i,ip]*R*unc[i]
                    sumξ += ξ[i,ip]
                end
            end
            Γ[i] = sumΓ
            dF1_dELL[i,i] = -1.0 + sumELL
            dF1_dRRP[i] = sumRRP
            dF1_dWWP[i] = sumWWP
            dF1_dTTP[i] = sumTTP
            dF1_dRR[i] = (ell[i] + sumξ*a[i])/R
            dF1_dWW[i] = sumξ*R*unc[i]*(1.0+1.0/ν)*sgrid_total[i]*η[i]
            dF1_dTT[i] = sumξ*R*unc[i]
        end
    end
    return dF1_dELLP, dF1_dRRP, dF1_dWWP, dF1_dTTP, dF1_dELL, dF1_dWW, dF1_dRR, dF1_dTT, unc, a
end

function kolmogorov_fwd_realbond(nx::Int, ns::Int,
                                 qfunction::Function, qp::Function,
                                 xgrid::Vector{Float64}, sgrid::Vector{Float64},
                                 ggrid::Vector{Float64}, xwts::Vector{Float64}, swts::Vector{Float64},
                                 unc::Vector{Float64}, a::Vector{Float64},
                                 R::Float64, ν::Float64, γ::Float64,
                                 ell::Vector{Float64}, μ::Vector{Float64}, η::Vector{Float64}, c::Vector{Float64})
    dF2_dMU   = zeros(nx*ns,nx*ns)
    dF2_dELL  = zeros(nx*ns,nx*ns)
    dF2_dRR = zeros(nx*ns)
    dF2_dWW = zeros(nx*ns)
    dF2_dTT = zeros(nx*ns)

    # KF jacobian terms
    for ixp = 1:nx
        for isp = 1:ns
            ip = ixp + nx*(isp-1)
            sumRR = 0.0
            sumWW = 0.0
            sumTT = 0.0
            for ix = 1:nx
                for is = 1:ns
                    i = ix + nx*(is-1)
                    mf = R*(sgrid[is]*η[i]+xgrid[ix] - c[i]) - xgrid[ixp]
                    dF2_dMU[ip,i]  = qfunction(mf)*ggrid[isp]*xwts[ix]*swts[is]*swts[isp]
                    dF2_dELL[ip,i] = qp(mf)*ggrid[isp]*xwts[ix]*swts[is]*swts[isp]*μ[ix]*R*unc[i]*((1/γ) + sgrid[is]*η[i]/(ν*c[i]))*ell[i]^(-(1/γ) - 1)
                    sumRR += qp(mf)*μ[ix]*ggrid[isp]*xwts[ix]*swts[is]*swts[isp]*unc[i]*a[i]/R
                    sumWW += qp(mf)*μ[ix]*ggrid[isp]*xwts[ix]*swts[is]*swts[isp]*unc[i]*R*(1+(1/ν))*sgrid[is]*η[i]
                    sumTT += qp(mf)*μ[ix]*ggrid[isp]*xwts[ix]*swts[is]*swts[isp]*unc[i]*R
                end
            end
            dF2_dRR[ip] = sumRR
            dF2_dWW[ip] = sumWW
            dF2_dTT[ip] = sumTT
        end
    end
    return dF2_dMU, dF2_dELL, dF2_dRR, dF2_dWW, dF2_dTT
end

function market_clearing_realbondmkup(nx::Int, ns::Int,
                                  xgrid::Vector{Float64}, sgrid::Vector{Float64},
                                  xwts::Vector{Float64}, swts::Vector{Float64},
                                  η::Vector{Float64}, c::Vector{Float64},
                                  μ::Vector{Float64}, ell::Vector{Float64},
                                  γ::Float64, ν::Float64,
                                  chipW::Vector{Float64}, chipR::Vector{Float64},
                                  chipX::Vector{Float64}, unc::Vector{Float64})
    # MKT clearing
    dF3_dMU  = zeros(nx*ns)
    dF3_dELL = zeros(nx*ns)
    dF3_dWW  = 0.0
    dF3_dRR  = 0.0
    dF3_dTT  = 0.0
    GDP      = 0.0

    for ix =1:nx
        for is=1:ns
            i = ix + nx*(is-1)
            dF3_dMU[i]  = xwts[ix]*swts[is]*(sgrid[is]*η[i] - c[i])
            dF3_dELL[i] = xwts[ix]*swts[is]*μ[ix]*(1+sgrid[is]*γ*η[i]/(ν*c[i]))*(1/γ)*unc[i]*ell[i]^(-(1/γ)-1)
            dF3_dWW += swts[is]*xwts[ix]*μ[i]*((sgrid[is]*η[i])/ν - (1+sgrid[is]*γ*η[i]/(ν*c[i]))*(1-unc[i])*chipW[i])
            dF3_dRR += -swts[is]*xwts[ix]*μ[i]*(1+sgrid[is]*γ*η[i]/(ν*c[i]))*(1-unc[i])*chipR[i]
            dF3_dTT += -swts[is]*xwts[ix]*μ[i]*(1+sgrid[is]*γ*η[i]/(ν*c[i]))*(1-unc[i])*chipX[i]
            GDP += swts[is]*xwts[ix]*sgrid[is]*η[i]*μ[i]
        end
    end

    return dF3_dMU, dF3_dELL, dF3_dWW, dF3_dRR, dF3_dTT, GDP
end

function normalize(m::RealBondMkup, JJ::Matrix{Float64})

    Qx, _, Qleft, Qright = compose_normalization_matrices(m)

    m <= Setting(:n_predetermined_variables, size(Qx, 1))

	Jac1 = Qleft*JJ*Qright

    return Jac1
end

function compose_normalization_matrices(m::RealBondMkup)
    nx::Int = get_setting(m, :nx)
    ns::Int = get_setting(m, :ns)

    # Create PPP matrix
    P1 = kron(eye(ns),ones(nx,1))
    P2 = kron(eye(ns), eye(nx)[:, 2:end])
    P  = hcat(P1, P2)

    Q,R = qr(P)
    Q = Array(Q)
    S         = Q[:, ns+1:end]'

    nxns = nx*ns

    Qleft     = cat(eye(nxns),S,[1],[1],[1],[1],[1],[1],[1],[1],[1], dims = [1 2])
    Qx        = cat(S,[1],[1],[1],[1], dims = [1 2])
    Qy        = cat(eye(nxns),[1],[1],[1],[1],[1], dims = [1 2])
    Qright    = cat(Qx',Qy',Qx',Qy', dims = [1 2])

    return Qx, Qy, Qleft, Qright
end
