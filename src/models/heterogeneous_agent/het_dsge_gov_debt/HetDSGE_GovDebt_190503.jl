using Distributions
#using PyPlot
using Roots, Random
using JLD2
Random.seed!(0)
save_data = true

function tauchen86(μ::AbstractFloat,ρ::AbstractFloat,σ::AbstractFloat,n::Int64,λ::AbstractFloat)
    #output is xgrid, xprob
    # x_t+1 = μ + ρ x_t + σ e_{t+1}, e_t+1 ∼ N(0,1)
    xhi = μ/(1-ρ) + λ*sqrt(σ^2/(1-ρ)^2)
    xlo = μ/(1-ρ) - λ*sqrt(σ^2/(1-ρ)^2)
    xgrid = zeros(n);
    xscale =(xhi-xlo)/(n-1)

    for i=1:n
        xgrid[i] = xlo + xscale*(i-1)
    end
    m=zeros(n-1);
    for i=1:n-1
        m[i] = (xgrid[i]+xgrid[i+1])/2
    end
    xprob = zeros(n,n) # xprob[i,j] = Pr(x_t+1=xgrid[j]|x_t=xgrid[i])
    for i=1:n # this is the state today
        normie = Normal(μ+ρ*xgrid[i],σ)
        normpdf(x) = pdf.(normie,x)
        normcdf(x) = cdf.(normie,x)
        for j=2:n-1
            xprob[i,j] = normcdf(m[j]) - normcdf(m[j-1])
        end
        xprob[i,1] = normcdf(m[1])
        xprob[i,n] = 1 - normcdf(m[n-1])
    end
    xprob = xprob./sum(xprob,2) # make sure the rows sum to 1
    return ( xgrid,xprob, xscale )
end

function trunc_lognpdf(x, UP, LO, mu, sig)
    # truncated log normal pdf
    logn = LogNormal(mu, sig)
    return (x-LO.>=-eps()).*(x-UP.<=eps()).*pdf.(logn,x)/(cdf.(logn,UP)-cdf.(logn,LO))
end


# set parameters
# PARAMETERS THAT AFFECT STEADY STATE:
# note, we have log utility, so no CRRA parameter
# and we have a quarterly model
r = 0.01                        # 4 percent steady state annual real interest rate
α = 0.3                         # capital share
H = 1.0                         # aggregate hours worked
δ  = 0.03                       # depreciation
sH_over_sL = 7.18182
pLH = 0.01                       # prob of going from low to high persistent skill
pHL = 0.0325                      # prob of going fom
γ = 0.#0.004                       # TFP growth
GoverY = 0.01#18                        # steady state govt spending/gdp
g = 1/(1-GoverY)                # steady state value of g_t, where G_t/Y_t = (1-(1/g_t))
η = 0.                         # borrowing constraint (normalized by TFP)
In    = 0.443993816237631       # normalizing constant for the mollifier
na = 300                         # cash on hand ditn grid points - set to smaller numbers for debugging
ns = 2# 5                         # skill distribution grid points
zlo = 0.145455                      # second income shock to mollify actual income
zhi = 2. - zlo                       # upper bound on this process
BoverY = 0.26

# AGGREGATE SHOCKS
ρB = 0.5 # persistence of discount factor shock
ρG = 0.5 # persistence of govt spending shock
ρZ = 0.5 # persistence of tfp growth shock
ρμ = 0.5 # persistence of investment shock
ρlamw = 0.5 # persistence of wage mkup shock
ρlamf = 0.5 # persistence of price mkup shock
ρmon = 0.5 # persistence of mon policy shock

# OTHER PARAMETERS THAT AFFECT DYNAMICS
spp  = 4.                        # second derivative of investment adjustment cost
lamw = 1.5                       # wage markup
ϕh   = 2.                        # inverse frisch elasticity
Φw   = 10.                       # rotemberg cost for wages
lamf = 1.5                       # price markup
Φp   = 1.                       # rotemberg cost for prices
ρR   = 0.75                      # persistence in taylor rule
ψπ   = 1.5                      # weight on inflation in taylor rule
ψy   = 0.5                       # weight on output growth in taylor rule
δb   = 1.                       # δb = 1 means balanced budget

κ_p  = .5
κ_w  = .5

# OPTIONS
TRUNCATE      = true # default: true
mindens       = 1e-8
RESCALE_XWTS  = true # default: true
SIMPLE_TAYLOR = false # default: false

# solve for steady state of aggregate scalar variables
# here I am assuming firms are subsidized in ss so that real mc = 1
Rk = r + δ                               # rental rate on capital
ω = (α^(α/(1-α)))*(1-α)*Rk^(-α/(1-α))    # real wage
kl = (α/(1-α))*(ω/Rk)*e^γ                # capital/labor ratio
k = kl*H                                 # capital
x = (1-(1-δ)*e^(-γ))*k                   # investment
y = (e^(-α*γ))*(k^α)*H^(1-α)             # gdp
bg = BoverY*y                            # govt debt
Tg = (e^(-γ) - 1/(1+r))*bg + (1-(1/g))*y # net lump sum taxes
T = Rk*k*e^(-γ) - x - Tg                 # net transfer to hhs

function persistent_skill_process(sH_over_sL::AbstractFloat, pLH::AbstractFloat,
                                  pHL::AbstractFloat, ns::Int)
    f1 = [[1-pLH pLH];[pHL 1-pHL]] # f1[i,j] is prob of going from i to j
    ss_skill_distr = [pHL/(pLH+pHL); pLH/(pLH+pHL)]
    slo = 1./(ss_skill_distr'*[1;sH_over_sL])
    sgrid = slo*[1;sH_over_sL]
    sscale = sgrid[2] - sgrid[1]
    swts     = (sscale/ns)*ones(ns) #quadrature weights
    f = f1./repmat(swts',ns,1)
    return (f, sgrid, swts)
end
(f, sgrid, swts) = persistent_skill_process(sH_over_sL, pLH, pHL, ns)

function cash_grid(sgrid::AbstractArray, ω::AbstractFloat, H::AbstractFloat, r::AbstractFloat,
                   η::AbstractFloat, γ::AbstractFloat,
                   T::AbstractFloat, zlo::AbstractFloat, na::Int)
    smin = minimum(sgrid)*zlo                           # lowest possible skill
    alo_ss = ω*smin*H - (1+r)*η*e^(-γ) + T + sgrid[1]*ω*H*0.05       # lowest possible cash on hand in ss

    alo = alo_ss                    # lower bound on cash on hand - could be < alo_ss
    ahi = max(alo*2, alo + 12.0)    # upper bound on cash on hand
    ascale = (ahi-alo)              # size of w grids

    # make grids
    agrid    = collect(linspace(alo,ahi,na)) #Evenly spaced grid
    awts = (ascale/na)*ones(na)          #quadrature weights
    return (agrid, awts)
end
(agrid, awts) = cash_grid(sgrid, ω, H, r, η, γ, T, zlo, na)

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

zgrid = linspace(zlo,zhi,na)
zwts = (zhi-zlo)/na
sumz=0.
for i=1:na
    sumz += mollifier(zgrid[i],zhi,zlo)*zwts
end

# experiment with a different q
ne = 100
σe = 0.01
(legrid, fe, sscale) = tauchen86(0.0,0.0,σe,ne,2.0)
egrid = exp.(legrid)
eprob = fe[1,:]

function convoluted_q(x::AbstractFloat, zhi::AbstractFloat, zlo::AbstractFloat, ne::Int, egrid::Vector{Float64}, eprob::Vector{Float64})
    sumne = 0.
    for i=1:ne
        sumne += mollifier(x - egrid[i], zhi, zlo)*eprob[i]
    end
    return sumne
end

qfunction(x) = mollifier(x, zhi, zlo)#/sumz # this ensures the mollifer sums to 1, at least on this particular grid
#qfunction(x) = convoluted_q(x, zhi, zlo, ne, egrid, eprob)

function dmollifier(x::AbstractFloat, ehi::AbstractFloat, elo::AbstractFloat)
    In = 0.443993816237631
    if x<ehi && x>elo
        temp = (-1.0 + 2.0*(x-elo)/(ehi-elo))
        out  = -(2*temp./((1 - temp.^2).^2)).*(2/(ehi-elo)).*mollifier(x, ehi, elo)
    else
        out = 0.0
    end
    return out
end


qp(z) = dmollifier(z,zhi,zlo)

Win = 1./(5.+0.02*(repeat(agrid,ns)-5))
Win = 2*ones(na*ns)/(agrid[end]+agrid[1])
aswts = kron(swts,awts)

function sspolicy(na::Int, ns::Int, β::AbstractFloat, R::AbstractFloat,
                     ω::AbstractFloat,
                     H::AbstractFloat,
                     η::AbstractFloat,
                     T::AbstractFloat,
                     γ::AbstractFloat,
                     qfunction::Function,
                     agrid::Vector{Float64},
                     sgrid::Vector{Float64},
                     aswts::Vector{Float64},
                     Win::Vector{Float64},
                     f::Array{Float64,2},damp::Float64 = 0.5, dist::Float64 = 1.,
                         tol::Float64 = 1e-4, maxit::Int64 = 500)
    n = na*ns
    c  = zeros(n)      # consumption
    bp = zeros(n)      # savings
    counter = 1
    Wout = copy(Win)
    while dist>tol && counter<maxit # for debugging
        # compute c(w) given guess for Win = β*R*E[u'(c_{t+1})]
        for iss in 1:ns
            for ia in 1:na
                c[na*(iss-1)+ia] =  min(1/Win[na*(iss-1)+ia],agrid[ia]+η)
            end
        end
        bp = repeat(agrid,ns) - c  # compute bp(w) given guess for Win
        Wout = parameterized_expectations(na,ns,β,R,ω,H,η,T,γ,qfunction,agrid,sgrid,aswts,c,bp,f)
        dist = maximum(abs.(Wout-Win))
        #println(dist)
        if mod(counter,25)==1
            println(dist)
            #plot(Wout)
        end
        Win = damp*Wout + (1.0-damp)*Win
        counter += 1
    end
    if counter==maxit
        warn("Euler iteration did not converge")
    end
    tr = kolmogorov_fwd(na,ns,ω,H,η,T,γ,qfunction,agrid,sgrid,bp,f)
    return (c,bp,Wout,tr)
end

function parameterized_expectations(na::Int,ns::Int, β::AbstractFloat, R::AbstractFloat, ω::AbstractFloat,
                                    H::AbstractFloat, η::AbstractFloat, T::AbstractFloat, γ::AbstractFloat,
                                    qfunction::Function,
                                    agrid::Vector{Float64}, sgrid::Vector{Float64},
                                    aswts::Vector{Float64}, c::Vector{Float64},
                                    bp::Vector{Float64}, f::Array{Float64,2})
    l_out = zeros(na*ns)
    for iss=1:ns
        for ia=1:na
            sumn = 0.0
            for isp=1:ns
                for iap=1:na
                    sumn += (aswts[na*(isp-1)+iap]/c[na*(isp-1)+iap])*qfunction((agrid[iap] - R*(e^(-γ))*bp[na*(iss-1)+ia] - T)/(ω*H*sgrid[isp]))*f[iss,isp]./sgrid[isp]
                end
            end
            l_out[na*(iss-1)+ia] = (β*R*(e^(-γ))/ω*H)*sumn
        end
    end
    return l_out
end

function kolmogorov_fwd(na::Int, ns::Int, ω::AbstractFloat,
                        H::AbstractFloat, η::AbstractFloat, T::AbstractFloat, γ::AbstractFloat,
                        qfunction::Function,agrid::Vector{Float64},
                        sgrid::Vector{Float64}, bp::Vector{Float64},f::Array{Float64,2})
    tr = zeros(na*ns,na*ns)
    for iss=1:ns
        for ia=1:na
            for isp=1:ns
                for iap=1:na
                    tr[na*(isp-1)+iap,na*(iss-1)+ia] = qfunction((agrid[iap] - R*(e^(-γ))*bp[na*(iss-1)+ia] - T)/(ω*H*sgrid[isp]))*f[iss,isp]./(ω*H*sgrid[isp])
                end
            end
        end
    end
    return tr
end

# just to test if code returns
R = 1+r
βlo = 0.5*(e^γ)/R # excess should be -ve
βhi = (e^γ)/R # excess should be +ve

β = (βlo+βhi)/2.0

function findss(na::Int, ns::Int, βlo::AbstractFloat,
                βhi::AbstractFloat, R::AbstractFloat,
                     ω::AbstractFloat,
                     H::AbstractFloat,
                     η::AbstractFloat,
                     T::AbstractFloat,
                     γ::AbstractFloat,
                     bg::AbstractFloat,
                     qfunction::Function,
                     agrid::Vector{Float64},
                     sgrid::Vector{Float64},
                     aswts::Vector{Float64},
                     Win::Vector{Float64},
                     f::Array{Float64,2},
                     excess::AbstractFloat = 5000., tol::AbstractFloat = 1e-4,
                     maxit::Int64 = 20)
    counter=1
    n       = ns*na
    c = zeros(n)
    bp = zeros(n)
    KF = zeros(n,n)
    m = zeros(n)
    report = zeros(maxit,2)
    β = (βlo+βhi)/2.0
    while abs(excess)>tol && counter<maxit # clearing markets
        β = (βlo+βhi)/2.0
        Win_guess = ones(n) # or Win
        (c, bp, Win, KF) = sspolicy(na, ns, β, R, ω, H, η, T, γ, qfunction, agrid, sgrid,
                                       aswts, Win_guess, f)

        LPMKF=aswts[1]*KF
        # find eigenvalue closest to 1
        (D,V) = eig(LPMKF)
        order_D = sortperm(abs.(D), rev = true)
        V = V[:,order_D]
        D = D[order_D]
        if abs(D[1]-1)>2e-1 # that's the tolerance we are allowing
            warn("your eigenvalue is too far from 1, something is wrong")
        end
        m = real(V[:,1]) #Pick the eigen vecor associated with the largest eigenvalue and moving it back to values
        m = m/(aswts'*m) #Scale of eigenvectors not determinate: rescale to integrate to exactly 1
        excess = (aswts'*(m.*bp))[1] - bg  #compute excess supply of savings, which is a fn of w
                # bisection
        println([βlo β βhi])
        println(excess)
        report[counter,1] = β # save the history of guesses and excess supply so we can check it's upward-sloping and monotonic
        report[counter,2] = excess
        if excess>0
            βhi=β
        elseif excess<0
            βlo = β
        end
        counter += 1
    end
    return (Win, c, m, β, report)
end

tic()
(ell, c, m, β, report) = findss(na, ns, βlo, βhi, R, ω, H, η, T, γ, bg, qfunction, agrid, sgrid,aswts, Win, f)
toc()

if save_data
    @show β
    file = JLD2.jldopen("steady_state.jld2", "w")
    # Aggregate variables
    write(file, "Rk", Rk)
    write(file, "omega", ω)
    write(file, "kl", kl)
    write(file, "k", k)
    write(file, "y", y)
    write(file, "T", T)

    # Functional/distributional variables
    write(file, "ell", ell)
    write(file, "c", c)
    write(file, "mu", m)
    write(file, "beta", β)
    close(file)
end


if TRUNCATE == true
    oldna = na
    na = maximum(find(m[1:na]+m[na+1:2*na].>mindens)) # used to be 1e-8
    m = m[[1:na;oldna+1:oldna+na]]
    ell = ell[[1:na;oldna+1:oldna+na]]
    # experiment with this:
    if RESCALE_XWTS==true
        agrid   = agrid[1:na] #Evenly spaced grid
        ahi = agrid[na]
        ascale = ahi-agrid[1]
        awts     = (ascale/na)*ones(na)          #quadrature weights
        aswts = kron(swts,awts)
    end
end

qp(z) = dmollifier(z,zhi,zlo)

n = na*ns
# we will always order things XP YP X Y. P denotes prime, i.e. t+1 variables
# convention is that capital letters generally refer to indices
# endogenous function-valued states:
KFP = 1:n           # combination of lagged ell function and lagged m function that predicts m
# endogenous scalar-valued states:
KKP   = n+1         # capital - we call the index KK to distinguish from the steady state object K
LRRP  = n+2         # lagged real interest rate
LIIP  = n+3        # lagged nominal rate
LYP   = n+4        # lagged gdp
LWP   = n+5        # lag real wages
LXP   = n+6        # lag investment
BGP   = n+7        # govt debt
# exogenous scalar-valued states:
BP    = n+8        # discount factor shock
GP    = n+9        # govt spending
ZP    = n+10       # tfp growth
MUP   = n+11        # investment shock
LAMWP = n+12        # wage markup
LAMFP = n+13        # price markup
MONP  = n+14        # monetary policy shock
# function-valued jumps
ELLP  = n+15:2*n+14 # ell function
# scalar-valued jumps
RRP   = 2*n+15        # real interest rate
IIP   = 2*n+16        # nominal interest rate
TTP   = 2*n+17        # transfers + dividends
WP    = 2*n+18        # real wage
HHP   = 2*n+19        # hours worked
PIP   = 2*n+20        # inflation
PIWP  = 2*n+21        # nominal wage inflation
LAMP  = 2*n+22        # average marginal utility
YP    = 2*n+23        # gdp
XP    = 2*n+24        # investment
MCP   = 2*n+25        # marginal cost - this is ζ in HetDSGE_kd.pdf
QP    = 2*n+26        # Tobin's qfunction
RKP   = 2*n+27        # return on capital
TGP   = 2*n+28        # lump sum tax
# now create all these indices also for current variables
nvars = 2*n+28
nscalars = 28 # number of eqs which output scalars
nyscalars = 14 # number of scalar jumps
nxscalars = nscalars - nyscalars # number of scalar states

KF = KFP + nvars
KK = KKP + nvars
LRR = LRRP + nvars
LII = LIIP + nvars
LY = LYP + nvars
LW = LWP + nvars
LX = LXP + nvars
BG = BGP + nvars
B = BP + nvars
G = GP + nvars
Z = ZP + nvars
MU = MUP + nvars
LAMW = LAMWP + nvars
LAMF = LAMFP + nvars
MON = MONP + nvars
ELL = ELLP + nvars
RR = RRP + nvars
II = IIP + nvars
TT = TTP + nvars
W = WP + nvars
HH = HHP + nvars
PI = PIP + nvars
PIW = PIWP + nvars
LAM = LAMP + nvars
Y = YP + nvars
X = XP + nvars
MC = MCP + nvars
Q = QP + nvars
RK = RKP + nvars
TG = TGP + nvars

# create objects needed for solve.jl
# we will order function blocks as follows:
# 1. all function blocks which output a function (first real eqs, then lags)
# 2. all function blocks which map functions to scalars
# 3. all scalar blocks involving endogenous vbls (first real eqs, then lags)
# 4. shock processes
funops = 1:2 # which operators output a function
# function blocks which output a function
F1  = 1:n # euler eqn
F2  = n+1:2*n # KF eqn
# function blocks which map functions to scalars
F5  = 2*n+1:2*n+1 # mkt clearing
F6  = 2*n+2:2*n+2 # lambda = average marginal utility
# scalar blocks involving endogenous variables
F7  = 2*n+3:2*n+3 # transfers
F8  = 2*n+4:2*n+4 # investment
F9  = 2*n+5:2*n+5 # tobin's q
F10 = 2*n+6:2*n+6 # capital accumulation
F11 = 2*n+7:2*n+7 # wage phillips curve
F12 = 2*n+8:2*n+8 # price phillips curve
F13 = 2*n+9:2*n+9 # marginal cost
F14 = 2*n+10:2*n+10 # gdp
F15 = 2*n+11:2*n+11 # optimal K/L ratio
F16 = 2*n+12:2*n+12 # taylor rule
F17 = 2*n+13:2*n+13 # fisher eqn
F18 = 2*n+14:2*n+14 # nominal wage inflation
F40 = 2*n+15:2*n+15 # fiscal rule
F41 = 2*n+16:2*n+16 # govt budget constraint
# lagged variables
F19 = 2*n+17:2*n+17 # LR
F20 = 2*n+18:2*n+18 # LI
F24 = 2*n+19:2*n+19 # LY
F31 = 2*n+20:2*n+20 # LW
F32 = 2*n+21:2*n+21 # LX
# shocks
F33 = 2*n+22:2*n+22 # discount factor B
F34 = 2*n+23:2*n+23 # govt spending G
F35 = 2*n+24:2*n+24 # tfp growth Z
F36 = 2*n+25:2*n+25 # investment MU
F37 = 2*n+26:2*n+26 # wage mkup LAMW
F38 = 2*n+27:2*n+27 # price mkup LAMF
F39 = 2*n+28:2*n+28 # monetary policy MON

# make the Jacobian
JJ = zeros(nvars, 2*nvars)

# first create some auxiliary terms for the eqns involving functions
# Euler equation
ee  = zeros(n,n) # ee[i,j] takes you from i to j
ξ   = zeros(n,n)
Ξ   = zeros(n,n)
unc = 1./ell .<= repeat(agrid,ns) + η
dF1_dELL = zeros(n,n)
dF1_dRZ  = zeros(n)
dF1_dELLP = zeros(n,n)
dF1_dWHP = zeros(n)
dF1_dTTP = zeros(n)
ellRHS = zeros(n)
for iss=1:ns
    for ia=1:na
        i  = na*(iss-1)+ia
        sumELL = 0.
        sumRZ  = 0.
        sumWH  = 0.
        sumTT  = 0.
        sum_ellRHS = 0.
        for isp=1:ns
            for iap=1:na
                ip = na*(isp-1)+iap
                ee[i,ip] = (agrid[iap] - R*(e^(-γ))*max(agrid[ia]-1/ell[i], -η) - T)/(ω*H*sgrid[isp])
                ξ[i,ip] = ((β*R*aswts[i]*e^(-γ))/(ω*H*sgrid[isp])^2)*max(ell[ip],1/(agrid[iap]+η))*qp(ee[i,ip])*f[iss,isp]
                Ξ[i,ip] = ((β*R*aswts[i]*e^(-γ))/(ω*H*sgrid[isp]))*max(ell[ip],1/(agrid[iap]+η))*qfunction(ee[i,ip])*f[iss,isp]
                sumELL += ξ[i,ip]*R*(e^(-γ))*unc[i]/ell[i]
                sumRZ  += ξ[i,ip]*R*(e^(-γ))*max(agrid[ia] - 1/ell[i],-η)
                dF1_dELLP[i,ip] = Ξ[i,ip]*unc[ip]
                sumWH  += Ξ[i,ip] + ξ[i,ip]*ee[i,ip]*(ω*H*sgrid[isp])
                sumTT  += ξ[i,ip]*T
                sum_ellRHS += Ξ[i,ip]
            end
        end
        ellRHS[i] = sum_ellRHS
        dF1_dELL[i,i] = -ell[i] - sumELL
        dF1_dRZ[i] = ellRHS[i] - sumRZ
        dF1_dRZ[i] = ell[i] - sumRZ
        dF1_dWHP[i] = -sumWH
        dF1_dTTP[i] = -sumTT
    end
end

# KF equation
bigΨ    = zeros(n,n) # this is also dF2_dM
smallψ = zeros(n,n)
dF2_dRZ = zeros(n)
dF2_dELL = zeros(n,n)
dF2_dM = zeros(n,n)
dF2_dWH = zeros(n)
dF2_dTT = zeros(n)

for isp=1:ns
    for iap=1:na
        ip  = na*(isp-1)+iap
        sumWH = 0.
        sumRZ = 0.
        sumTT = 0.
        for iss=1:ns
            for ia=1:na
                i = na*(iss-1)+ia
                bigΨ[ip,i] = aswts[i]*m[i]*qfunction(ee[i,ip])*f[iss,isp]/(ω*H*sgrid[isp])
                smallψ[ip,i] = aswts[i]*m[i]*qp(ee[i,ip])*f[iss,isp]/((ω*H*sgrid[isp])^2)
                sumWH += bigΨ[ip,i] + smallψ[ip,i]*ee[i,ip]*(ω*H*sgrid[isp])
                sumRZ += smallψ[ip,i]*(R*e^(-γ))*max(agrid[ia] - 1/ell[i],-η)
                dF2_dELL[ip,i] = smallψ[ip,i]*(R*e^(-γ))*(unc[i]/ell[i])
                sumTT += smallψ[ip,i]*T
                dF2_dM[ip,i] = aswts[i]*qfunction(ee[i,ip])*f[iss,isp]/(ω*H*sgrid[isp]) # note, now we linearize
            end
        end
        dF2_dWH[ip] = -sumWH
        dF2_dRZ[ip] = -sumRZ
        dF2_dTT[ip] = -sumTT
    end
end

#mkt clearing, lambda function
c = min.(1./ell,repeat(agrid,ns)+η)
lam = (aswts.*m)'*(1./c) # average marginal utility which the union uses to set wages
ϕ = lam*ω/(H^ϕh) # now that we know lam in steady state, choose disutility to target hours H

# Euler equation
JJ[F1,ELLP] = dF1_dELLP
JJ[F1,ZP]   = -dF1_dRZ
JJ[F1,WP]   = dF1_dWHP
JJ[F1,HHP]  = dF1_dWHP
JJ[F1,TTP]  = dF1_dTTP
JJ[F1,B]    = ell
JJ[F1,ELL]  = dF1_dELL
JJ[F1,RR]   = dF1_dRZ

# KF eqn
JJ[F2,KFP]  = -eye(n)
JJ[F2,KF]   = dF2_dM
JJ[F2,ELL]  = dF2_dELL
JJ[F2,RR]    = dF2_dM*dF2_dRZ
JJ[F2,Z]    = -dF2_dM*dF2_dRZ
JJ[F2,W]    = dF2_dM*dF2_dWH
JJ[F2,HH]   = dF2_dM*dF2_dWH
JJ[F2,TT]   = dF2_dM*dF2_dTT

# mkt clearing
JJ[F5,Y]   = y/g
JJ[F5,G]   = -y/g
JJ[F5,X]   = -x
JJ[F5,ELL] = (m.*unc.*aswts.*c)'
JJ[F5,KF]   = -(aswts.*c)' # note, now we linearize
JJ[F5,RR]   = -(aswts.*c)'*dF2_dRZ
JJ[F5,Z]   = (aswts.*c)'*dF2_dRZ
JJ[F5,W]   = -(aswts.*c)'*dF2_dWH
JJ[F5,HH]   = -(aswts.*c)'*dF2_dWH
JJ[F5,TT]   = -(aswts.*c)'*dF2_dTT

# lambda = average marginal utility
JJ[F6,LAM] = lam
JJ[F6,KF]   = -(aswts./c)' # note, now we linearize
JJ[F6,RR]   = -(aswts./c)'*dF2_dRZ
JJ[F6,Z]   = (aswts./c)'*dF2_dRZ
JJ[F6,W]   = -(aswts./c)'*dF2_dWH
JJ[F6,HH]   = -(aswts./c)'*dF2_dWH
JJ[F6,TT]   = -(aswts./c)'*dF2_dTT
JJ[F6,ELL] = -(aswts.*unc.*m./c)'

# transfer
JJ[F7,TT]  = T
JJ[F7,RK]  = -Rk*k
JJ[F7,KK]  = -Rk*k
JJ[F7,Z]   = Rk*k
JJ[F7,X]   = x
JJ[F7,MC]  = y
JJ[F7,TG]  = Tg

# investment
JJ[F8,Q]  = 1.
JJ[F8,MU] = 1.
JJ[F8,XP] = spp*(e^(3*γ))/R
JJ[F8,ZP] = spp*(e^(3*γ))/R
JJ[F8,X]  = -spp*(e^(3*γ))/R - spp*e^(2*γ)
JJ[F8,LX] = spp*e^(2*γ)
JJ[F8,Z]  = -spp*e^(2*γ)

# tobin's q
JJ[F9,RR]  = R
JJ[F9,Q]   = R
JJ[F9,RKP] = -Rk
JJ[F9,QP]  = -(1-δ)

# capital accumulation
JJ[F10,KKP] = 1.
JJ[F10,KK]  = -(1-δ)
JJ[F10,Z]   = (1-δ)
JJ[F10,MU]  = -x/k
JJ[F10,X]   = -x/k

# wage phillips curve
JJ[F11,PIW]  = -1.
JJ[F11,LAMW] = 1. #(ϕ*H^ϕh)/Φw
JJ[F11,HH]   = κ_w*ϕh #(ϕ*(H^ϕh)*(1+lamw)/lamw*Φw)*ϕh
JJ[F11,W]    = -κ_w #-(ϕ*(H^ϕh)*(1+lamw)/lamw*Φw)
JJ[F11,PIWP]  = β
JJ[F11,LAM]  = -κ_w #-(ϕ*(H^ϕh)*(1+lamw)/lamw*Φw)

# price phillips curve
JJ[F12,PI]   = -1.
JJ[F12,MC]   = κ_p #(1+lamf)/(lamf*Φp)
JJ[F12,LAMF] = 1. #1/Φp
JJ[F12,PIP]  = 1/R

# marginal cost
JJ[F13,MC] = 1.
JJ[F13,W]  = -(1-α)
JJ[F13,RK] = -α

# gdp
JJ[F14,Y]  = 1.
JJ[F14,Z]  = α
JJ[F14,KK] = -α
JJ[F14,HH] = -(1-α)

# optimal k/l ratio
JJ[F15,RK] = 1.
JJ[F15,W]  = -1.
JJ[F15,HH] = -1.
JJ[F15,KK] = 1.
JJ[F15,Z]  = -1.

# taylor rule
JJ[F16,II]   = -1.
JJ[F16,LII]  = ρR
JJ[F16,PI]   = (1-ρR)*ψπ
JJ[F16,Y]    = (1-ρR)*ψy
JJ[F16,LY]  = -(1-ρR)*ψy
JJ[F16,Z]    = (1-ρR)*ψy
JJ[F16,MON] = 1.

# fisher eqn
JJ[F17,RR]  = 1.
JJ[F17,PIP] = 1.
JJ[F17,II]  = -1.

# wage inflation
JJ[F18,PIW] = 1.
JJ[F18,PI]  = -1.
JJ[F18,Z]   = -1.
JJ[F18,W]   = -1.
JJ[F18,LW]  = 1.

# fiscal rule
JJ[F40,TG]  = -Tg
JJ[F40,RR]   = δb*bg/R
JJ[F40,BG]  = δb*bg*e^(-γ)
JJ[F40,Z]   = -δb*bg*e^(-γ)
JJ[F40,Y]   = δb*(1-(1/g))*y
JJ[F40,G]   = δb*y/g


# govt budget constraint
JJ[F41,BGP] = -(bg/R)
JJ[F41,RR]   = bg/R
JJ[F41,BG]  = bg*e^(-γ)
JJ[F41,Z]   = -bg*e^(-γ)
JJ[F41,Y]   = (1-(1/g))*y
JJ[F41,G]   = y/g
JJ[F41,TG]  = -Tg

# update lagged variables
JJ[F19,LRRP] = 1.
JJ[F19,RR]   = -1.

JJ[F20,LIIP] = 1.
JJ[F20,II]   = -1.

JJ[F24,LYP] = 1.
JJ[F24,Y]   = -1.

JJ[F31,LWP] = 1.
JJ[F31,W]   = -1.

JJ[F32,LXP] = 1.
JJ[F32,X]   = -1.

# discount factor shock
JJ[F33,BP] = 1.
JJ[F33,B]  = -ρB

# g/y shock
JJ[F34,GP] = 1.
JJ[F34,G]  = -ρG

# tfp growth shock
JJ[F35,ZP] = 1.
JJ[F35,Z]  = -ρZ

# investment shock
JJ[F36,MUP] = 1.
JJ[F36,MU]  = -ρμ

# wage mkup shock
JJ[F37,LAMWP] = 1.
JJ[F37,LAMW]  = -ρlamw

# price mkup shock
JJ[F38,LAMFP] = 1.
JJ[F38,LAMF]  = -ρlamf

# monetary policy shock
JJ[F39,MONP] = 1.
JJ[F39,MON]  = -ρmon

if save_data
    JLD2.jldopen("jacobian.jld2", "w") do file
        file["JJ"] = JJ
        file["KFP"] = KFP
        file["KKP"]   = KKP
        file["LRRP"]  = LRRP
        file["LIIP"]  = LIIP
        file["LYP"]   = LYP
        file["LWP"]   = LWP
        file["LXP"]   = LXP
        file["BP"]    = BP
        file["GP"]    = GP
        file["ZP"]    = ZP
        file["MUP"]   = MUP
        file["LAMWP"] = LAMWP
        file["LAMFP"] = LAMFP
        file["MONP"]  = MONP
        file["ELLP"]  = ELLP
        file["RRP"]   = RRP
        file["IIP"]   = IIP
        file["TTP"]   = TTP
        file["WP"]    = WP
        file["HHP"]   = HHP
        file["PIP"]   = PIP
        file["PIWP"]  = PIWP
        file["LAMP"]  = LAMP
        file["YP"]   = YP
        file["XP"]    = XP
        file["MCP"]   = MCP
        file["QP"]    = QP
        file["RKP"]   = RKP
	file["KF"]    = KF
        file["KK"] = KK
        file["LRR"] = LRR
        file["LII"] = LII
        file["LY"] = LY
        file["LW"] = LW
        file["LX"] = LX
        file["B"] = B
        file["G"] = G
        file["Z"] = Z
        file["MU"] = MU
        file["LAMW"] = LAMW
        file["LAMF"] = LAMF
        file["MON"] = MON
        file["ELL"] = ELL
        file["RR"] = RR
        file["II"] = II
        file["TT"] = TT
        file["W"] = W
        file["HH"] = HH
        file["PI"] = PI
        file["PIW"] = PIW
        file["LAM"] = LAM
        file["Y"] = Y
        file["X"] = X
        file["MC"] = MC
        file["Q"] = Q
        file["RK"] = RK
        file["TGP"] = TGP
        file["BGP"] = BGP
        file["TG"] = TG
        file["BG"] = BG
        file["F1"]  = F1
        file["F2"]  = F2
        file["F5"]  = F5
        file["F6"]  = F6
        file["F7"]  = F7
        file["F8"]  = F8
        file["F9"]  = F9
        file["F10"] = F10
        file["F11"] = F11
        file["F12"] = F12
        file["F13"] = F13
        file["F14"] = F14
        file["F15"] = F15
        file["F16"] = F16
        file["F17"] = F17
        file["F18"] = F18
        file["F19"] = F19
        file["F20"] = F20
        file["F24"] = F24
        file["F31"] = F31
        file["F32"] = F32
        file["F33"] = F33
        file["F34"] = F34
        file["F35"] = F35
        file["F36"] = F36
        file["F37"] = F37
        file["F38"] = F38
        file["F39"] = F39
        file["F40"] = F40
        file["F41"] = F41
    end
end

# create matrix used to normalize distributions
P1 = kron(eye(ns),ones(na,1))
Ptemp = eye(na)
Ptemp = Ptemp[:,2:end]
P2 = kron(eye(ns),Ptemp)
P = [P1 P2]

(QQQ,Rjunk)=qr(P)

S         = QQQ[:,ns+1:end]'; #
Qleft     = cat([1 2],eye(n),S,eye(nscalars))
Qx        = cat([1 2],S,eye(nxscalars))
Qy        = cat([1 2],eye(n),eye(nyscalars))

include("solve.jl")

tic()
gx2, hx2, gx, hx =  solve(JJ, Qleft, Qx, Qy)
toc()

if save_data
    JLD2.jldopen("solve.jld2", "w") do file
        file["gx"] = gx
        file["hx"] = hx
        file["gx2"] = gx2
        file["hx2"] = hx2
    end
end
