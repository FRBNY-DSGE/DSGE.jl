function jacobian(m::BondLabor)

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
    chipW = (1+ν)./( ν./( aborrow + χss - xgrid_total) + γ./χss  )
    chipR = (1/(R*R))*(  ν*abar./(aborrow + χss - xgrid_total) )./( ν./( aborrow + χss - xgrid_total) + γ./χss  )

    # we will always order things XP YP X Y
    # convention is that capital letters generally refer to indices
    MUP  = 1:nx*ns
    ZP   = nx*ns+1
    ELLP = nx*ns+2:2*nx*ns+1
    RP   = 2*nx*ns+2

    MU   = 2*nx*ns+3:3*nx*ns+2
    Z    = 3*nx*ns+3
    ELL  = 3*nx*ns+4:4*nx*ns+3
    R    = 4*nx*ns+4

    # create objects needed for solve.jl
    F1 = 1:nx*ns # euler eqn
    F2 = nx*ns+1:2*nx*ns # KF
    F3 = 2*nx*ns+1:2*nx*ns+1 # mkt ckr
    F4 = 2*nx*ns+2:2*nx*ns+2 # z

    # EE
    ξ   = zeros(nx*ns,nx*ns)
    Γ   = zeros(nx*ns,1)
    ellRHS = zeros(nx*ns,1)
    Ξ   = zeros(nx*ns,nx*ns)
    LEE = zeros(nx*ns,nx*ns)
    WEE = zeros(nx*ns,1)
    REE = zeros(nx*ns,1)

    # Mollifier convenience functions
    qfunction(x) = mollifier(x,ehi,elo)
    qp(x) = dmollifier(x,ehi,elo)

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
        end
    end

    ELLEE = (β/γ)*R*R*Γ.*(1+ (γ/ν)*sgrid_total.*η./c).*(l.^(-1-1/γ) ).*((l.^(-1/γ)).<=χss)-1

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
    WMKT = 0.0
    RMKT = 0.0

    for ix in 1:nx
        for is in 1:ns
            i = ix + nx*(is-1)
            RMKT += -( swts[is]*μ[i]*γ*η[i]/(ν*c[i]) + μ[i] )*swts[is]*xwts[ix]*chipR[i]*((l[i]^(-1/γ))>χss[i])
            WMKT += (1+1/ν)*η[i]*sgrid[is]*μ[i]*swts[is]*xwts[ix] - (1+ sgrid[is]*γ*η[i]/(ν*c[i]))*μ[i]*swts[is]*xwts[ix]*chipW[i]*((l[i]^(-1/γ))>χss[i])
        end
    end

    # Make the Jacobian
    JJ = zeros(2*nx*ns+2,4*nx*ns+4)

    # Euler equation
    JJ[F1,ZP]  = -γ*β*R*WEE

    JJ[F1,ELLP] = β*R*LEE

    JJ[F1,RP] = -γ*β*R*REE

    JJ[F1,Z]  = β*R*R*Γ.*(sgrid_total.*η*(1+1/ν) - (1+ (γ/ν)*sgrid_total.*η./c).*chipW.*((l.^(-1/γ)).>χss)  )

    JJ[F1,ELL] = diagm(ELLEE[:])

    JJ[F1,R] = ellRHS/R + β*R*(sgrid_total.*η + xgrid_total - c).*Γ - β*R*R*Γ.*(1+ (γ/ν)*sgrid_total.*η./c ).*chipR.*((l.^(-1/γ)).>χss)

    # KF
    JJ[F2,MUP] = eye(nx*ns)

    JJ[F2,MU]  = MKF

    JJ[F2,Z]   = WKF

    JJ[F2,ELL] = LKF

    JJ[F2,R]   = RKF

    # mkt ckr
    JJ[F3,MU]  = (η'.*sgrid_total' - c').*weights_total'

    JJ[F3,Z]   = WMKT'

    JJ[F3,ELL] = (1/γ)*μ'.*(sgrid_total'*(γ/ν).*(η./c)' + 1).*weights_total'.*( l.^(-1-1/γ) )'.*((l.^(-1/γ)).<=χss)'

    JJ[F3,R]   = RMKT'

    # TFP
    JJ[F4,ZP]  = -1

    JJ[F4,Z]   = ρ_z

    return JJ
end
