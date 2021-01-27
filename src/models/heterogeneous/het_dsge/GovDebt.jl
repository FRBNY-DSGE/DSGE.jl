using Distributions
using PyPlot
using Roots
using JLD

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
    xprob = xprob ./ sum(xprob,2) # make sure the rows sum to 1
    return ( xgrid,xprob, xscale )
end

function trunc_lognpdf(x, UP, LO, mu, sig)
    # truncated log normal pdf
    logn = LogNormal(mu, sig)
    return (x-LO.>=-eps()) .* (x-UP.<=eps()).*pdf.(logn,x)/(cdf.(logn,UP)-cdf.(logn,LO))
end


# set parameters
# PARAMETERS THAT AFFECT STEADY STATE:
# note, we have log utility, so no CRRA parameter
# and we have a quarterly model
r = 0.01                        # 4 percent steady state annual real interest rate
α = 0.3                         # capital share
H = 1.0                         # aggregate hours worked
δ  = 0.03                       # depreciation
sH_over_sL = 1.2 / 0.8
pLH = 0.1                       # prob of going from low to high persistent skill
pHL = 0.1                      # prob of going fom
γ = 0.#0.004                       # TFP growth
GoverY = 0.01#18                        # steady state govt spending/gdp
g = 1/(1-GoverY)                # steady state value of g_t, where G_t/Y_t = (1-(1/g_t))
η = 0.                         # borrowing constraint (normalized by TFP)
In    = 0.443993816237631       # normalizing constant for the mollifier
na = 300                         # cash on hand ditn grid points - set to smaller numbers for debugging
ns = 2# 5                         # skill distribution grid points
zlo =  1. / 3.                      # second income shock to mollify actual income
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

# OPTIONS
TRUNCATE      = true # default: true
mindens       = 1e-8
RESCALE_XWTS  = true # default: true
SIMPLE_TAYLOR = false # default: false

# solve for steady state of aggregate scalar variables
# here I am assuming firms are subsidized in ss so that real mc = 1
Rk = r + δ                               # rental rate on capital
ω = (α^(α/(1-α)))*(1-α)*Rk^(-α/(1-α))    # real wage
kl = (α/(1-α))*(ω/Rk)*ℯ^γ                # capital/labor ratio
k = kl*H                                 # capital
x = (1-(1-δ)*ℯ^(-γ))*k                   # investment
y = (ℯ^(-α*γ))*(k^α)*H^(1-α)             # gdp
bg = BoverY*y                            # govt debt
T = Rk*k*ℯ^(-γ) - x - (1-(1/g))*y - ((1+r)*ℯ^(-γ)-1)*bg        # net transfer to hhs

function persistent_skill_process(sH_over_sL::AbstractFloat, pLH::AbstractFloat, pHL::AbstractFloat, ns::Int)
    f1 = [[1-pLH pLH];[pHL 1-pHL]] # f1[i,j] is prob of going from i to j
    ss_skill_distr = [pHL/(pLH+pHL); pLH/(pLH+pHL)]
    slo = 1. / (ss_skill_distr'*[1;sH_over_sL])
    sgrid = slo*[1;sH_over_sL]
    sscale = sgrid[2] - sgrid[1]
    swts     = (sscale/ns)*ones(ns) #quadrature weights
    f = f1 ./ repmat(swts',ns,1)
    return (f, sgrid, swts)
end
(f, sgrid, swts) = persistent_skill_process(sH_over_sL, pLH, pHL, ns)

function cash_grid(sgrid::AbstractArray, ω::AbstractFloat, H::AbstractFloat, r::AbstractFloat, η::AbstractFloat, γ::AbstractFloat,
                   T::AbstractFloat, zlo::AbstractFloat, na::Int)
    smin = minimum(sgrid)*zlo                           # lowest possible skill
    alo_ss = ω*smin*H - (1+r)*η*ℯ^(-γ) + T + sgrid[1]*ω*H*0.05       # lowest possible cash on hand in ss

    alo = alo_ss                    # lower bound on cash on hand - could be < alo_ss
    ahi = max(alo*2,alo+5.)                      # upper bound on cash on hand
    ascale = (ahi-alo)              # size of w grids

    # make grids
    agrid    = collect(linspace(alo,ahi,na)) #Evenly spaced grid
    awts = (ascale/na)*ones(na)          #quadrature weights
    return (agrid, awts)
end
smin = minimum(sgrid)*zlo                           # lowest possible skill
alo_ss = ω*smin*H - (1+r)*η*ℯ^(-γ) + T + sgrid[1]*ω*H*0.05       # lowest possible cash on hand in ss

alo = alo_ss                    # lower bound on cash on hand - could be < alo_ss
ahi = max(alo*2,alo+5.)                      # upper bound on cash on hand
ascale = (ahi-alo)              # size of w grids

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
        out  = -(2*temp ./ ((1 - temp.^2).^2)) .* (2/(ehi-elo)).*mollifier(x, ehi, elo)
    else
        out = 0.0
    end
    return out
end


qp(z) = dmollifier(z,zhi,zlo)

Win = 1 ./ (5. + 0.02*(repeat(agrid,ns)-5))
Win = 2*ones(na*ns)/(ahi+alo)
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
                    sumn += (aswts[na*(isp-1)+iap]/c[na*(isp-1)+iap])*qfunction((agrid[iap] - R*(ℯ^(-γ))*bp[na*(iss-1)+ia] - T)/(ω*H*sgrid[isp]))*f[iss,isp] ./ sgrid[isp]
                end
            end
            l_out[na*(iss-1)+ia] = (β*R*(ℯ^(-γ))/ω*H)*sumn
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
                    tr[na*(isp-1)+iap,na*(iss-1)+ia] = qfunction((agrid[iap] - R*(ℯ^(-γ))*bp[na*(iss-1)+ia] - T)/(ω*H*sgrid[isp]))*f[iss,isp] ./ (ω*H*sgrid[isp])
                end
            end
        end
    end
    return tr
end

# just to test if code returns
R = 1+r
βlo = 0.5*(ℯ^γ)/R # excess should be -ve
βhi = (ℯ^γ)/R # excess should be +ve

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
    #Win0=copy(Win)
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
        excess = (aswts'*(m .* bp))[1] - bg  #compute excess supply of savings, which is a fn of w
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
