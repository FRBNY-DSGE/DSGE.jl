using JLD
using Distributions
using Roots
include("klein_solve.jl")

test_output = true

function tauchen86(μ::AbstractFloat,σ::AbstractFloat,n::Int64,λ::AbstractFloat)
    #output is xgrid, xprob
    xhi = μ + λ*σ;
    xlo = μ - λ*σ;
    xgrid = zeros(n);
    xscale =(xhi-xlo)/(n-1)

    for i=1:n
        xgrid[i] = xlo + xscale*(i-1);
    end
    m=zeros(n-1);
    for i=1:n-1
        m[i] = (xgrid[i]+xgrid[i+1])/2;
    end
    xprob = zeros(n);
    normie = Normal(μ,σ)
    normpdf(x) = pdf.(normie,x)
    normcdf(x) = cdf.(normie,x)
    for j=2:n-1
        xprob[j] = normcdf(m[j]) - normcdf(m[j-1]);
    end
    xprob[1] = normcdf(m[1]);
    xprob[n] = 1 - normcdf(m[n-1]);
    return ( xgrid,xprob, xscale )
end

# set parameters
# Economic Parameters
R        = 1.04                    # 4 percent steady state real interst rate
γ        = 1.0                     # CRRA Parameter
ν        = 1.0                     # Inverse Frisch Elasticity of Labor supply
abar     = -0.5                    # Borrowing floor
ρz       = 0.95                    # persistence of agg tfp
mus      = 0.0
sigs     = 0.5                     # sigma of log normal in income

# Settings
In       = 0.443993816237631       # normalizing constant for the mollifier
nx       = 50                      # cash on hand ditn grid points
ns       = 2                      # skill distribution grid points
n        = nx*ns
# slo      = 0.5                     # lower bound on skill grid
# shi      = 1.5                     # upper bound on skill
# sscale   = shi-slo                 # size of our s grid
elo      = 0.0                     # stochastic consumption commitments
ehi      = 1.0
xlo      = abar - ehi              # lower bound on cash on hand
xhi      = 4.0                    # upper bound on cash on hand
xscale   = xhi-xlo                 # size of w grid

# make grids
xgrid    = collect(linspace(xlo,xhi,nx)) # Evenly spaced grid
xwts     = (xscale/nx)*ones(nx)          # quadrature weights

if ns == 1
    sgrid = ones(1)
    swts = ones(1)
    g=ones(1)
else
    (lsgrid,sprob,sscale) = tauchen86(mus,sigs,ns,3.0)
    swts  = (sscale/ns)*ones(ns)    # quadrature weights
    sgrid = exp.(lsgrid)            # exp(sgrid)
    g     = sprob./swts
end

# make kron grid, for stacking nx*ns as opposed to dealing with two dimensions
SGRIDB   = kron(sgrid,ones(nx))
XGRIDB   = kron(ones(ns),xgrid)
BIGWTS   = kron(swts,xwts)

# Normalizing constants
# Note: xscale/nx == xwts[1], and sscale/ns == swts[1]
MX       = sqrt(xscale/nx)*eye(nx)       # matrix that maps from set of values of a functions at a point to the set of normalized scaling function coefficients
MS       = sqrt(sscale/ns)*eye(ns)
MXinv    = sqrt(nx/(xscale))*eye(nx)     # inverse of MX
MSinv    = sqrt(ns/(sscale))*eye(ns)     # inverse of MS
MSMX = kron(MS,MX)
MSMXinv = inv(MSMX)

function mollifier(z::AbstractFloat,ehi::AbstractFloat,elo::AbstractFloat)
    # mollifier function
    In = 0.443993816237631
    if z<ehi && z>elo
        temp = -1.0 + 2.0*(z-elo)/(ehi-elo)
        out = (2.0/(ehi-elo))*exp(-1.0/(1.0-temp*temp))/In
    else
        out = 0.0
    end
    return out
end

function dmollifier(x::AbstractFloat, ehi::AbstractFloat, elo::AbstractFloat)
    In = 0.443993816237631
    temp = (-1.0 + 2.0*(x-elo)/(ehi-elo))
    dy  = -(2*temp./((1 - temp.^2).^2)).*(2/(ehi-elo)).*mollifier(x, ehi, elo)
end

qfunction(x) = mollifier(x,ehi,elo)
qp(x) = dmollifier(x,ehi,elo)

# guess for the W function W(w) = beta R E u'(c_{t+1})
Win = ones(nx,ns)

# compute constrained consumption and labor supply once and for all
χss = zeros(nx*ns)
# ηss = zeros(nx*ns)
aborrow  = abar/R

laborzero(s::AbstractFloat,x::AbstractFloat,c::AbstractFloat,ν::AbstractFloat,γ::AbstractFloat,aborrow::AbstractFloat) = s^(1+ν) - ((aborrow+c-x)^ν)*c^γ
chifun(s::AbstractFloat,x::AbstractFloat,ν::AbstractFloat,γ::AbstractFloat,aborrow::AbstractFloat) = fzero(c ->laborzero(s,x,c,ν,γ,aborrow), 0.0,15.0)
# etafun(s::AbstractFloat,x::AbstractFloat,ν::AbstractFloat,γ::AbstractFloat,aborrow::AbstractFloat) = chifun(s,x,ν,γ,aborrow)^(-(γ/ν))*s^(1.0/ν)

for iss=1:ns
    for ix=1:nx
        χss[ix + nx*(iss-1)] = chifun(sgrid[iss],xgrid[ix],ν,γ,aborrow)
        # ηss[ix + nx*(iss-1)] = etafun(sgrid[iss],xgrid[ix],ν,γ,aborrow)
    end
end

chipW = (1+ν)./( ν./( aborrow + χss - XGRIDB) + γ./χss  )
chipR = (1/(R*R))*(  ν*abar./(aborrow + χss - XGRIDB) )./( ν./( aborrow + χss - XGRIDB) + γ./χss  )

# initial guesses
βguess = 0.8

winconst= βguess*R*sum((qfunction.(abar - XGRIDB).*χss.^(-γ)).*(kron((swts.*g),xwts)))
# println(winconst)
Win = winconst*ones(nx*ns)

function sspolicy(β::AbstractFloat,
                  R::AbstractFloat,
                  γ::AbstractFloat,
                  sigs::AbstractFloat,
                  ehi::AbstractFloat,
                  elo::AbstractFloat,
                  XGRIDB::Vector{Float64},
                  SGRIDB::Vector{Float64},
                  xwts::Vector{Float64},
                  swts::Vector{Float64},
                  Win::Vector{Float64},
                  g::Vector{Float64},
                  χss::Vector{Float64})
    # outputs: [c,ap,Wout,tr]

    nx      = length(xwts)
    ns      = length(swts)
    n       = ns*nx
    c       = zeros(n)      # consumption
    h       = zeros(n)      # hours
    ap      = zeros(n)      # savings
    damp    = 0.05          # how much we update W guess
    counter = 1
    dist    = 1.0          # distance between Win and Wout on decision rules
    tol     = 1e-4          # acceptable error
    Wout    = copy(Win)    # initialize Wout
    tr      = zeros(n,n)   # transition matrix mapping todays dist of w to w'
    qxg = zeros(n,n)          # auxilary variable to compute LHS of euler
    maxit = 500
    while dist>tol && counter<maxit
        # compute c(w) given guess for Win = β*R*E[u'(c_{t+1})]
        c  = min.(Win.^(-1.0/γ),χss)
        h  = (SGRIDB.^(1.0/ν)).*(c.^(-γ/ν))
        ap = R*(XGRIDB+SGRIDB.*h - c)
        for ix=1:nx
            for iss=1:ns
                sumn=0.0
                for ixp=1:nx
                    for isp = 1:ns
                       sumn += qfunction(ap[ix+nx*(iss-1)] - XGRIDB[nx*(isp-1)+ixp])*g[isp]*xwts[ixp]*swts[isp]*(c[nx*(isp-1)+ixp]^(-γ) )
                    end
                end
                Wout[ix+nx*(iss-1)] = β*R*sumn
            end
        end
        dist = maximum(abs.(Wout-Win))
        # if counter%1==0
        #     println(dist)
        # end
        Win = damp*Wout + (1.0-damp)*Win
        counter += 1
    end
    if counter==maxit
        warn("Euler iteration did not converge")
    end
    # c  = min.(Win.^(-1.0/γ),χss)
    # h  = (SGRIDB.^(1.0/ν)).*(c.^(-γ/ν))
    # ap = R*(XGRIDB+SGRIDB.*h - c)

    for ix=1:nx
        for ixp=1:nx
            for iss = 1:ns
                for isp=1:ns
                    tr[ixp+nx*(isp-1),ix+nx*(iss-1)] += qfunction(ap[ix+nx*(iss-1)] - XGRIDB[nx*(isp-1)+ixp])*g[isp]
                end
            end
        end
    end
    return (c,h,ap,Wout,tr)
end

βlo = 0.4/R # excess should be -ve
βhi = 1/R # excess is +
function findss(βlo::AbstractFloat,
                βhi::AbstractFloat,
                R::AbstractFloat,
                γ::AbstractFloat,
                sigs::AbstractFloat,
                ehi::AbstractFloat,
                elo::AbstractFloat,
                XGRIDB::Vector,
                SGRIDB::Vector,
                xwts::Vector,
                swts::Vector,
                Win::Vector{Float64},
                g::Vector{Float64},
                χss::Vector{Float64})
    excess = 5000.0 # excess supply of savings
    tol = 1e-5
    count=1
    maxit=10
    nx      = length(xwts)
    ns      = length(swts)
    n       = ns*nx
    c = zeros(n)
    h = zeros(n)
    ap = zeros(n)
    KF = zeros(n,n)
    μ = zeros(n)
    β = (βlo+βhi)/2.0
    xswts = kron(swts,xwts)
    while abs(excess)>tol && count<maxit # clearing markets
        # println(excess)
        β = (βlo+βhi)/2.0
        (c, h, ap, Win, KF) = sspolicy(β, R, γ, sigs,
                                       ehi, elo, XGRIDB, SGRIDB,
                                       xwts, swts, Win, g, χss)
        # sspolicy returns the consumption and hours decision rules, a prime
        # decision rule, Win = updated β R u'(c_{t+1}),
        # KF is the Kolmogorov foward operator
        # and is the map which moves you from cash on hand distribution
        # today to cash on had dist tomorrow

        LPMKF = MSMX*KF*MSMX'

        temp = xwts[1]*swts[1]*KF
        @assert LPMKF ≈ temp

        # find eigenvalue closest to 1
        (Dee,Vee) = eig(LPMKF)
        if abs(Dee[1]-1)>2e-1 # that's the tolerance we are allowing
            warn("your eigenvalue is ", Dee[1], " which is too far from 1, something is wrong")
        end

        μ = real(Vee[:,1]) # Pick the eigen vecor associated with the largest
        # eigenvalue and moving it back to values

        μ = μ/(xswts'*μ) # Scale of eigenvectors not determinate: rescale to integrate to exactly 1
        # change all following stuff to work with R not K

        excess = xswts'*(μ.*ap)  # compute excess supply of savings, which is a fn of w

        # bisection
        if excess>0
            βhi=β
        elseif excess<0
            βlo = β
        end
        count += 1
    end
    return (Win, c, h, μ, β)
end

tic()
(ell, c, η, μ, β) = findss(βlo, βhi, R, γ, sigs,
                           ehi, elo, XGRIDB, SGRIDB,
                           xwts, swts, Win, g,χss)
toc()

test_output && include("test/steady_state.jl")

# Jacobian stuff
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

# Create auxiliary variables
# c = min.(ell.^(-1.0/γ),χss) # Consumption Decision Rule
# η = (SGRIDB.^(1.0/ν)).*(c.^(-γ/ν))

# EE
ξ   = zeros(nx*ns,nx*ns)
Γ   = zeros(nx*ns,1)
ellRHS = zeros(nx*ns,1)
Ξ   = zeros(nx*ns,nx*ns)
LEE = zeros(nx*ns,nx*ns)
WEE = zeros(nx*ns,1)
REE = zeros(nx*ns,1)

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

tic()
gx2, hx2, gx, hx =  klein_solve(JJ, Qleft, Qx, Qy)
toc()

test_output && include("test/solve.jl")
