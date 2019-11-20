function jacobian(m::BondLabor)

    # Load in endogenous state and eq cond indices
    endo = augment_model_states(m.endogenous_states_unnormalized,
                                n_model_states_unnormalized(m))
    eq   = m.equilibrium_conditions

    # Load in parameters, steady-state parameters, and grids
    γ    = m[:γ]
    abar = m[:abar]
    R    = m[:R]
    ν    = m[:ν]
    ρ_z  = m[:ρ_z].value

    l    = m[:lstar].value
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

    # Jacobian stuff

    # Create auxiliary variables
    # c = min.(l.^(-1.0/γ),χss) # Consumption Decision Rule
    # η = (sgrid_total.^(1.0/ν)).*(c.^(-γ/ν))
    aborrow  = abar/R
    chipW = (1+ν)./( ν./( aborrow .+ χss - xgrid_total) + γ./χss  )
    chipR = (1/(R*R))*(  ν*abar./(aborrow .+ χss - xgrid_total) )./( ν./( aborrow .+ χss - xgrid_total) + γ./χss  )

    # EE
    ξ   = zeros(nx*ns,nx*ns)
    Γ   = zeros(nx*ns,1)
    ellRHS = zeros(nx*ns,1)
    Ξ   = zeros(nx*ns,nx*ns)
    LEE = zeros(nx*ns,nx*ns)
    WEE = zeros(nx*ns,1)
    REE = zeros(nx*ns,1)
    ELLEE = zeros(nx*ns)
    ELWEE = zeros(nx*ns)
    ELREE = zeros(nx*ns)

    # Mollifier convenience functions
    qfunction(x) = mollifier_bondlabor(x,ehi,elo)
    qp(x) = dmollifier_bondlabor(x,ehi,elo)

    # EE jacobian terms
    for ix in 1:nx
        for is in 1:ns
            i = ix + nx*(is-1)
            sumΓ   = 0.0
            sumWEE = 0.0
            sumREE = 0.0
            sumELLRHS = 0.0
            for ixp in 1:nx
                for isp in 1:ns
                    ip = ixp + nx*(isp-1)
                    mf = R*(sgrid[is]*η[i]+xgrid[ix] - c[i]) - xgrid[ixp]
                    ξ[i,ip] = qp(mf)*ggrid[isp]
                    Ξ[i,ip] = (c[ip]^(-γ-1))*qfunction(mf)*ggrid[isp]*swts[isp]*xwts[ixp]
                    LEE[i,ip] = Ξ[i,ip]*(l[ip]^(-1- 1/γ))*((l[ip]^(-1/γ))<=χss[ip])
                    sumΓ+= (c[ip]^(-γ))*ξ[i,ip]*swts[isp]*xwts[ixp]
                    sumWEE += Ξ[i,ip]*chipW[ip]*((l[ip]^(-1/γ))>χss[ip])
                    sumREE += Ξ[i,ip]*chipR[ip]*((l[ip]^(-1/γ))>χss[ip])
                    sumELLRHS += β*R*(c[ip]^(-γ))*qfunction(mf)*ggrid[isp]*swts[isp]*xwts[ixp]
                end
            end
            Γ[i] = sumΓ
            WEE[i] = sumWEE
            REE[i] = sumREE
            ellRHS[i] = sumELLRHS

            ELWEE[i] = sumΓ * (sgrid_total[i]*η[i]*(1+1/ν) - (1+(γ/ν)*sgrid_total[i]*η[i]/c[i])*chipW[i]*((l[i]^(-1/γ))>χss[i]))
            ELREE[i] = β*R*(sgrid_total[i] * η[i] + xgrid_total[i] - c[i])*sumΓ - β*R^2*sumΓ*(1 + (γ/ν)*sgrid_total[i]*η[i]/c[i]) * chipR[i] * ((l[i]^(-1/γ))>χss[i])
            ELLEE[i] = (β/γ)*R^2*sumΓ * (1 + (γ/ν)*sgrid_total[i]*η[i]/c[i]) * (l[i]^(-1-1/γ)) * ((l[i]^(-1/γ)) <= χss[i])-1
        end
    end

    # KF
    LKF = zeros(nx*ns,nx*ns)
    WKF = zeros(nx*ns,1)
    RKF = zeros(nx*ns,1)
    MKF = zeros(nx*ns,nx*ns)

    # KF jacobian terms
    for ixp in 1:nx
        for isp in 1:ns
            ip = ixp + nx*(isp-1)
            sumWKF1 = 0.0
            sumWKF2 = 0.0
            sumRKF = 0.0
            for ix in 1:nx
                for is in 1:ns
                    i = ix + nx*(is-1)
                    mf = R*(sgrid[is]*η[i]+xgrid[ix] - c[i]) - xgrid[ixp]
                    MKF[ip,i] = qfunction(mf)*ggrid[isp]*swts[is]*xwts[ix]
                    LKF[ip,i] = R*ξ[i,ip]*μ[i]*( sgrid[is]*γ*η[i]/(ν*c[i]) + 1 )*swts[is]*xwts[ix]*(l[i]^(-1/γ-1))*((l[i]^(-1/γ))<=χss[i])/γ
                    sumWKF1 += R*(1+ 1/ν)*ξ[i,ip]*μ[i]*sgrid[is]*η[i]*swts[is]*xwts[ix]
                    sumWKF2 += -R*ξ[i,ip]*μ[i]*( sgrid[is]*γ*η[i]/(ν*c[i]) + 1 )*swts[is]*xwts[ix]*chipW[i]*((l[i]^(-1/γ))>χss[i])
                    sumRKF  += ξ[i,ip]*μ[i]*( sgrid[is]*η[i] + xgrid[ix] - c[i] - R*(sgrid[is]*γ*η[i]/(ν*c[i]) + 1 )*chipR[i]*((l[i]^(-1/γ))>χss[i]) )*swts[is]*xwts[ix]
                end
            end
            WKF[ip] = sumWKF1  + sumWKF2
            RKF[ip] = sumRKF
        end
    end

    # MKT CKR
    μMKT = zeros(nx*ns)
    LMKT = zeros(nx*ns)
    RMKT = 0.0
    WMKT = 0.0

    for ix in 1:nx
        for is in 1:ns
            i = ix + nx*(is-1)
            μMKT[i] = (η[i] * sgrid_total[i] - c[i]) * weights_total[i]
            LMKT[i] = (1/γ)*μ[i] * (sgrid_total[i]*(γ/ν)*(η[i]/c[i]) + 1) * weights_total[i] * (l[i]^(-1-1/γ))*((l[i]^(-1/γ))<=χss[i])
            RMKT += -( sgrid[is]*μ[i]*γ*η[i]/(ν*c[i]) + μ[i] )*swts[is]*xwts[ix]*chipR[i]*((l[i]^(-1/γ))>χss[i])
            WMKT += (1+1/ν)*η[i]*sgrid[is]*μ[i]*swts[is]*xwts[ix] - (1+ sgrid[is]*γ*η[i]/(ν*c[i]))*μ[i]*swts[is]*xwts[ix]*chipW[i]*((l[i]^(-1/γ))>χss[i])
        end
    end

    # Make the Jacobian
    JJ = zeros(2*nx*ns+2,4*nx*ns+4)

    # Euler equation
    JJ[eq[:eq_euler], endo[:z′_t]]  = -γ*β*R*WEE

    JJ[eq[:eq_euler], endo[:l′_t]] = β*R*LEE

    JJ[eq[:eq_euler], endo[:R′_t]] = -γ*β*R*REE

    JJ[eq[:eq_euler], endo[:z_t]]  = β*R^2*ELWEE

    JJ[eq[:eq_euler], endo[:l_t]] = Diagonal(ELLEE)

    JJ[eq[:eq_euler], endo[:R_t]] = ellRHS/R + ELREE

    # KF
    JJ[eq[:eq_kolmogorov_fwd], endo[:μ′_t]] = eye(nx*ns)

    JJ[eq[:eq_kolmogorov_fwd], endo[:μ_t]]  = MKF

    JJ[eq[:eq_kolmogorov_fwd], endo[:z_t]]   = WKF

    JJ[eq[:eq_kolmogorov_fwd], endo[:l_t]] = LKF

    JJ[eq[:eq_kolmogorov_fwd], endo[:R_t]]   = RKF

    # mkt ckr
    JJ[eq[:eq_market_clearing], endo[:μ_t]]  =  μMKT

    JJ[eq[:eq_market_clearing], endo[:z_t]]  .= WMKT

    JJ[eq[:eq_market_clearing], endo[:l_t]]  =  LMKT

    JJ[eq[:eq_market_clearing], endo[:R_t]]  .= RMKT

    # TFP
    JJ[eq[:eq_TFP], endo[:z′_t]]  .= -1

    JJ[eq[:eq_TFP], endo[:z_t]]   .= ρ_z

    if !m.testing && get_setting(m, :normalize_distr_variables)
        JJ = normalize(m, JJ)
    end

    return JJ
end

function normalize(m::BondLabor, JJ::Matrix{Float64})

    Qx, _, Qleft, Qright = compose_normalization_matrices(m)

    m <= Setting(:n_predetermined_variables, size(Qx, 1))

	Jac1 = Qleft*JJ*Qright

    return Jac1
end

function compose_normalization_matrices(m::BondLabor)
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
    Qleft     = cat(S,[1],eye(nx*ns),[1], dims = [1 2])
    Qx        = cat(S,[1], dims = [1 2])
    Qy        = cat(eye(nx*ns),[1], dims = [1 2])
    Qright    = cat(Qx',Qy',Qx',Qy', dims = [1 2])

    return Qx, Qy, Qleft, Qright
end
