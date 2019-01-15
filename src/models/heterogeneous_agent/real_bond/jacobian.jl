function jacobian(m::RealBond)

    # Load in endogenous state and eq cond indices
    endo = augment_model_states(m.endogenous_states_unnormalized,
                                n_model_states_unnormalized(m))
    eq   = m.equilibrium_conditions

    # Load in parameters, steady-state parameters, and grids
    γ     = m[:γ]
    abar  = m[:abar]
    R     = m[:R]
    ν     = m[:ν]
    ρz    = m[:ρ_z].value
    ρmon  = m[:ρmon].value
    κ     = m[:κ].value
    phipi = m[:phipi].value
    aborrow = abar/R

    ell    = m[:lstar].value
    c    = m[:cstar].value
    μ    = m[:μstar].value
    η    = m[:ηstar].value
    β    = m[:βstar].value
    χss  = m[:χstar].value

    xgrid = m.grids[:xgrid].points
    xwts  = m.grids[:xgrid].weights
    sgrid = m.grids[:sgrid].points
    swts  = m.grids[:sgrid].weights
    ggrid = m.grids[:ggrid]
    xgrid_total = m.grids[:xgrid_total]
    sgrid_total = m.grids[:sgrid_total]
    weights_total = m.grids[:weights_total]

    elo = get_setting(m, :elo)
    ehi = get_setting(m, :ehi)

    nx = get_setting(m, :nx)
    ns = get_setting(m, :ns)

    # Mollifiers
    qfunction(x) = mollifier_realbond(x,ehi,elo)
    qp(x) = dmollifier_realbond(x,ehi,elo)

    # Jacobian stuff
    chipW = (1+ν)./(ν./(aborrow + χss - xgrid_total) + γ./χss)
    chipR = (1/(R*R))*(ν*abar./(aborrow + χss - xgrid_total))./(ν./(aborrow + χss - xgrid_total) + γ./χss)
    chipX = (ν*χss)./(ν*χss + γ*(aborrow + χss - xgrid_total))

    # we will always order things XP YP X Y
    # convention is that capital letters generally refer to indices
    MUP  = 1:nx*ns
    ZP   = nx*ns+1
    MONP = nx*ns+2 # monetary policy shock
    ELLP = nx*ns+3:2*nx*ns+2
    RRP  = 2*nx*ns+3
    IIP  = 2*nx*ns+4 # nominal rates
    WWP  = 2*nx*ns+5 # wages
    PIP  = 2*nx*ns+6 # inflation
    TTP  = 2*nx*ns+7 # transfers
    MU   = 2*nx*ns+8:3*nx*ns+7
    Z    = 3*nx*ns+8
    MON  = 3*nx*ns+9
    ELL  = 3*nx*ns+10:4*nx*ns+9
    RR   = 4*nx*ns+10
    II   = 4*nx*ns+11
    WW   = 4*nx*ns+12
    PI   = 4*nx*ns+13
    TT   = 4*nx*ns+14

    # create objects needed for solve.jl
    funops = 1:2 # which operators output a function
    F1 = 1:nx*ns             # euler eqn
    F2 = nx*ns+1:2*nx*ns     # KF
    F3 = 2*nx*ns+1:2*nx*ns+1 # mkt ckr
    F4 = 2*nx*ns+2:2*nx*ns+2 # z
    F5 = 2*nx*ns+3:2*nx*ns+3 # phillips
    F6 = 2*nx*ns+4:2*nx*ns+4 # taylor
    F7 = 2*nx*ns+5:2*nx*ns+5 # fisher
    F8 = 2*nx*ns+6:2*nx*ns+6 # transfers
    F9 = 2*nx*ns+7:2*nx*ns+7 # monetary policy

    # Create auxiliary variables
    # c = min.(ell.^(-1.0/γ),χss) # Consumption Decision Rule
    # η = (sgrid_total.^(1.0/ν)).*(c.^(-γ/ν))

    # EE
    ξ   = zeros(nx*ns,nx*ns)
    Γ   = zeros(nx*ns,1)
    ellRHS = zeros(nx*ns,1)
    Ξ   = zeros(nx*ns,nx*ns)
    dF1_dELLP = zeros(nx*ns,nx*ns)
    dF1_dRRP  = zeros(nx*ns,1)
    dF1_dWWP  = zeros(nx*ns,1)
    dF1_dTTP  = zeros(nx*ns,1)
    dF1_dELL  = zeros(nx*ns,nx*ns)
    dF1_dWW   = zeros(nx*ns,1)
    dF1_dRR   = zeros(nx*ns,1)
    dF1_dTT   = zeros(nx*ns,1)
    unc = zeros(nx*ns,1)
    a   = zeros(nx*ns,1)

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

    #KF
    dF2_dMU   = zeros(nx*ns,nx*ns)
    dF2_dELL = zeros(nx*ns,nx*ns)
    dF2_dRR = zeros(nx*ns,1)
    dF2_dWW = zeros(nx*ns,1)
    dF2_dTT = zeros(nx*ns,1)

    #KF jacobian terms
    for ixp=1:nx
        for isp =1:ns
            ip = ixp + nx*(isp-1)
            sumRR = 0.0
            sumWW = 0.0
            sumTT = 0.0
            for ix=1:nx
                for is = 1:ns
                    i = ix + nx*(is-1)
                    mf = R*(sgrid[is]*η[i]+xgrid[ix] - c[i]) - xgrid[ixp]
                    dF2_dMU[ip,i]  = qfunction(mf)*ggrid[isp]*xwts[ix]*swts[is]*swts[isp]
                    dF2_dELL[ip,i] = qp(mf)*ggrid[isp]*xwts[ix]*swts[is]*swts[isp]*μ[ix]*ggrid[isp]*R*unc[i]*((1/γ) + sgrid[is]*η[i]/(ν*c[i]))
                    sumRR +=qp(mf)*μ[ix]*ggrid[isp]*xwts[ix]*swts[is]*swts[isp]*unc[i]*a[i]/R
                    sumWW +=qp(mf)*μ[ix]*ggrid[isp]*xwts[ix]*swts[is]*swts[isp]*unc[i]*R*(1+(1/ν))*sgrid[is]*η[i]
                    sumTT +=qp(mf)*μ[ix]*ggrid[isp]*xwts[ix]*swts[is]*swts[isp]*unc[i]*R
                end
            end
            dF2_dRR[ip] = sumRR
            dF2_dWW[ip] = sumWW
            dF2_dTT[ip] = sumTT
        end
    end

    # MKT clearing
    dF3_dMU  = zeros(1,nx*ns)
    dF3_dELL = zeros(1,nx*ns)
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

    # Make the Jacobian
    JJ = zeros(2*nx*ns+7,4*nx*ns+14)

    # Euler equation
    JJ[F1,ELLP] = dF1_dELLP
    JJ[F1,RRP]  = dF1_dRRP
    JJ[F1,WWP]  = dF1_dWWP
    JJ[F1,TTP]  = dF1_dTTP
    JJ[F1,ELL]  = dF1_dELL
    JJ[F1,RR]   = dF1_dRR
    JJ[F1,WW]   = dF1_dWW
    JJ[F1,TT]   = dF1_dTT

    # KF
    JJ[F2,MUP] = -eye(nx*ns)
    JJ[F2,MU]  = dF2_dMU
    JJ[F2,ELL] = dF2_dELL
    JJ[F2,WW]  = dF2_dWW
    JJ[F2,RR]  = dF2_dRR
    JJ[F2,TT]  = dF2_dTT

    # mkt ckr
    JJ[F3,MU]  = dF3_dMU
    JJ[F3,Z]   = GDP
    JJ[F3,ELL] = dF3_dELL
    JJ[F3,WW]  = dF3_dWW
    JJ[F3,RR]  = dF3_dRR
    JJ[F3,TT]  = dF3_dTT

    # TFP
    JJ[F4,ZP]  = 1.0
    JJ[F4,Z]   = -ρz

    # Phillips
    JJ[F5,PIP] = -1.0/R
    JJ[F5,Z]   = κ
    JJ[F5,WW]  = -κ
    JJ[F5,PI]  = 1.0

    # taylor
    JJ[F6,II]  = 1.0
    JJ[F6,PI]  = -phipi*R
    JJ[F6,MON] = -1.0

    # fisher
    JJ[F7,II]  = 1.0
    JJ[F7,PIP] = -R
    JJ[F7,RR]  = -1.0

    # transfers
    JJ[F8,TT]  = 1.0
    JJ[F8,Z]   = -GDP
    JJ[F8,WW]  = GDP

    # MP
    JJ[F9,MONP] = 1.0
    JJ[F9,MON]  = -ρmon

    if !m.testing && get_setting(m, :normalize_distr_variables)
        JJ = normalize(m, JJ)
    end

    return JJ
end

function normalize(m::RealBond, JJ::Matrix{Float64})

    Qx, _, Qleft, Qright = compose_normalization_matrices(m)

    m <= Setting(:n_predetermined_variables, size(Qx, 1))

	Jac1 = Qleft*JJ*Qright

    return Jac1
end

function compose_normalization_matrices(m::RealBond)
    nx = get_setting(m, :nx)
    ns = get_setting(m, :ns)

    # Create PPP matrix
    P1 = kron(eye(ns),ones(nx,1))
    Ptemp = eye(nx)
    Ptemp = Ptemp[:,2:end]
    P2 = kron(eye(ns),Ptemp)
    P = [P1 P2]

    (Q,R) = qr(P)

    S         = Q[:,ns+1:end]'
    Qleft     = cat([1 2],eye(nx*ns),S,[1],[1],[1],[1],[1],[1],[1])
    Qx        = cat([1 2],S,[1],[1])
    Qy        = cat([1 2],eye(nx*ns),[1],[1],[1],[1],[1])
	Qright    = cat([1,2],Qx',Qy',Qx',Qy')

    return Qx, Qy, Qleft, Qright
end
