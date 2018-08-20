using DSGE
using BenchmarkTools
using Distributions
using JLD
using Roots

include("klein_solve.jl")

test_output = true

m = BondLabor()

# Steady-state computation
@btime steadystate!(m)

test_output && include("test/steady_state.jl")

throw("stop!")

# Jacobian stuff

# Create auxiliary variables
# c = min.(ell.^(-1.0/γ),χss) # Consumption Decision Rule
# η = (SGRIDB.^(1.0/ν)).*(c.^(-γ/ν))
aborrow  = abar/R
chipW = (1+ν)./( ν./( aborrow + χss - XGRIDB) + γ./χss  )
chipR = (1/(R*R))*(  ν*abar./(aborrow + χss - XGRIDB) )./( ν./( aborrow + χss - XGRIDB) + γ./χss  )

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
funops = 1:2 # which operators output a function
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
for ix=1:nx
    for is =1:ns
        i = ix + nx*(is-1)
        sumΓ   = 0.0
        sumWEE = 0.0
        sumREE = 0.0
        sumELLRHS = 0.0
        for ixp=1:nx
            for isp = 1:ns
                ip = ixp + nx*(isp-1)
                mf = R*(sgrid[is]*η[i]+xgrid[ix] - c[i]) - xgrid[ixp]
                ξ[i,ip] = qp(mf)*g[isp]
                Ξ[i,ip] = (c[ip]^(-γ-1))*qfunction(mf)*g[isp]*swts[isp]*xwts[ixp]
                LEE[i,ip] = Ξ[i,ip]*(ell[ip]^(-1- 1/γ))*((ell[ip]^(-1/γ))<=χss[ip])
                sumΓ+= (c[ip]^(-γ))*ξ[i,ip]*swts[isp]*xwts[ixp]
                sumWEE += Ξ[i,ip]*chipW[ip]*((ell[ip]^(-1/γ))>χss[ip])
                sumREE += Ξ[i,ip]*chipR[ip]*((ell[ip]^(-1/γ))>χss[ip])
                sumELLRHS += β*R*(c[ip]^(-γ))*qfunction(mf)*g[isp]*swts[isp]*xwts[ixp]
            end
        end
        Γ[i] = sumΓ
        WEE[i] = sumWEE
        REE[i] = sumREE
        ellRHS[i] = sumELLRHS
    end
end

ELLEE = (β/γ)*R*R*Γ.*(1+ (γ/ν)*SGRIDB.*η./c).*(ell.^(-1-1/γ) ).*((ell.^(-1/γ)).<=χss)-1

# KF
LKF = zeros(nx*ns,nx*ns)
WKF = zeros(nx*ns,1)
RKF = zeros(nx*ns,1)
MKF = zeros(nx*ns,nx*ns)

# KF jacobian terms
for ixp=1:nx
    for isp =1:ns
        ip = ixp + nx*(isp-1)
        sumWKF1 = 0.0
        sumWKF2 = 0.0
        sumRKF = 0.0
        for ix=1:nx
            for is = 1:ns
                i = ix + nx*(is-1)
                mf = R*(sgrid[is]*η[i]+xgrid[ix] - c[i]) - xgrid[ixp]
                MKF[ip,i] = qfunction(mf)*g[isp]*swts[is]*xwts[ix]
                LKF[ip,i] = R*ξ[i,ip]*μ[i]*( sgrid[is]*γ*η[i]/(ν*c[i]) + 1 )*swts[is]*xwts[ix]*(ell[i]^(-1/γ-1))*((ell[i]^(-1/γ))<=χss[i])/γ
                sumWKF1 += R*(1+ 1/ν)*ξ[i,ip]*μ[i]*sgrid[is]*η[i]*swts[is]*xwts[ix]
                sumWKF2 += -R*ξ[i,ip]*μ[i]*( sgrid[is]*γ*η[i]/(ν*c[i]) + 1 )*swts[is]*xwts[ix]*chipW[i]*((ell[i]^(-1/γ))>χss[i])
                sumRKF  += ξ[i,ip]*μ[i]*( sgrid[is]*η[i] + xgrid[ix] - c[i] - R*(sgrid[is]*γ*η[i]/(ν*c[i]) + 1 )*chipR[i]*((ell[i]^(-1/γ))>χss[i]) )*swts[is]*xwts[ix]
            end
        end
        WKF[ip] = sumWKF1  + sumWKF2
        RKF[ip] = sumRKF
    end
end

# MKT CKR
WMKT = 0.0
RMKT = 0.0

for ix =1:nx
    for is=1:ns
        i = ix + nx*(is-1)
        RMKT += -( swts[is]*μ[i]*γ*η[i]/(ν*c[i]) + μ[i] )*swts[is]*xwts[ix]*chipR[i]*((ell[i]^(-1/γ))>χss[i])
        WMKT += (1+1/ν)*η[i]*sgrid[is]*μ[i]*swts[is]*xwts[ix] - (1+ sgrid[is]*γ*η[i]/(ν*c[i]))*μ[i]*swts[is]*xwts[ix]*chipW[i]*((ell[i]^(-1/γ))>χss[i])
    end
end

# Make the Jacobian
JJ = zeros(2*nx*ns+2,4*nx*ns+4)

# Euler equation
JJ[F1,ZP]  = -γ*β*R*WEE

JJ[F1,ELLP] = β*R*LEE

JJ[F1,RP] = -γ*β*R*REE

JJ[F1,Z]  = β*R*R*Γ.*(SGRIDB.*η*(1+1/ν) - (1+ (γ/ν)*SGRIDB.*η./c).*chipW.*((ell.^(-1/γ)).>χss)  )

JJ[F1,ELL] = diagm(ELLEE[:])

JJ[F1,R] = ellRHS/R + β*R*(SGRIDB.*η + XGRIDB - c).*Γ - β*R*R*Γ.*(1+ (γ/ν)*SGRIDB.*η./c ).*chipR.*((ell.^(-1/γ)).>χss)

# KF
JJ[F2,MUP] = eye(nx*ns)

JJ[F2,MU]  = MKF

JJ[F2,Z]   = WKF

JJ[F2,ELL] = LKF

JJ[F2,R]   = RKF

# mkt ckr
JJ[F3,MU]  = (η'.*SGRIDB' - c').*BIGWTS'

JJ[F3,Z]   = WMKT'

JJ[F3,ELL] = (1/γ)*μ'.*(SGRIDB'*(γ/ν).*(η./c)' + 1).*BIGWTS'.*( ell.^(-1-1/γ) )'.*((ell.^(-1/γ)).<=χss)'

JJ[F3,R]   = RMKT'

# TFP
JJ[F4,ZP]  = -1

JJ[F4,Z]   = ρz

test_output && include("test/jacobian.jl")

# Create PPP matrix
P1 = kron(eye(ns),ones(nx,1))
Ptemp = eye(nx)
Ptemp = Ptemp[:,2:end]
P2 = kron(eye(ns),Ptemp)
P = [P1 P2]

(Q,R)=qr(P)

S         = Q[:,ns+1:end]'; #
Qleft     = cat([1 2],S,[1],eye(nx*ns),[1])
Qx        = cat([1 2],S,[1])
Qy        = cat([1 2],eye(nx*ns),[1])

# Solve
tic()
gx2, hx2, gx, hx =  klein_solve(JJ, Qleft, Qx, Qy)
toc()

test_output && include("test/solve.jl")
