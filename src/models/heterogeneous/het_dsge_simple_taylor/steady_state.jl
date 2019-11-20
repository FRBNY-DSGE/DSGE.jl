function steadystate!(m::HetDSGESimpleTaylor;
                      βlo::Float64 = 0.5*exp(m[:γ])/(1 + m[:r]),
                      βhi::Float64 = exp(m[:γ])/(1 + m[:r]),
                      excess::Float64 = 5000.,
                      tol::Float64 = 1e-4,
                      maxit::Int64 = 20)

    # Load settings
    nx = get_setting(m, :nx)
    ns = get_setting(m, :ns)
    xlo = get_setting(m, :xlo)
    xhi = get_setting(m, :xhi)
    zlo = get_setting(m, :zlo)
    zhi = get_setting(m, :zhi)

    # Load parameters
    R = 1 + m[:r]
    H = m[:H].value
    η = m[:η].value
    γ = m[:γ].value
    ω = m[:ωstar].value
    T = m[:Tstar].value

    # Load grids
    xgrid = m.grids[:xgrid].points
    sgrid = m.grids[:sgrid].points
    xswts = m.grids[:weights_total]
    f     = m.grids[:fgrid]

    # Construct qfunction arguments
    zgrid = range(zlo, stop = zhi, length = nx)
    zwts  = (zhi-zlo)/nx
    sumz  = 0.
    for i in 1:nx
        sumz += mollifier_hetdsge(zgrid[i],zhi,zlo)*zwts
    end

    if get_setting(m, :use_last_βstar) && !isnan(m[:βstar].value)
        βlo = βhi = m[:βstar]
    end

    counter = 1
    n = ns*nx
    c = zeros(n)
    bp = zeros(n)
    KF = zeros(n,n)
    μ = zeros(n)

    # Initial guess
    β   = 1.0
    Win = Vector{Float64}(undef, n)
    while abs(excess) > tol && counter < maxit # clearing markets
        β = (βlo+βhi)/2.0
        Win_guess = ones(n)

        c, bp, Win, KF = policy_hetdsge(nx, ns, β, R, ω, H, η, T, γ, zhi, zlo, sumz, xgrid, sgrid,
                                          xswts, Win_guess, f)

        LPMKF = xswts[1]*KF

        # find eigenvalue closest to 1
        (D,V) = (eigen(LPMKF)...,)
        order_D = sortperm(abs.(D), rev = true)
        V = V[:,order_D]
        D = D[order_D]
        if abs(D[1]-1)>2e-1 # that's the tolerance we are allowing
            @warn "your eigenvalue is too far from 1, something is wrong"
        end

        μ = real(V[:,1]) # Pick the eigen vector associated with the largest eigenvalue and moving it back to values
        μ = μ/(xswts'*μ) # Scale of eigenvectors not determinate: rescale to integrate to exactly 1
        excess = (xswts'*(μ.*bp))[1]  # compute excess supply of savings, which is a fn of w
                                      # bisection
        if excess > 0
            βhi = β
        elseif excess < 0
            βlo = β
        end
        counter += 1
    end

    m[:lstar] = Win
    m[:cstar] = c
    m[:μstar] = μ
    m[:βstar] = β

    nothing
end

function policy_hetdsge(nx::Int, ns::Int, β::AbstractFloat, R::AbstractFloat,
                        ω::AbstractFloat, H::AbstractFloat,
                        η::AbstractFloat, T::AbstractFloat,
                        γ::AbstractFloat, zhi::AbstractFloat,
                        zlo::AbstractFloat, sumz::AbstractFloat,
                        xgrid::Vector{Float64}, sgrid::Vector{Float64},
                        xswts::Vector{Float64}, Win::Vector{Float64},
                        f::Array{Float64,2}, damp::Float64 = 0.5, dist::Float64 = 1.,
                        tol::Float64 = 1e-4, maxit::Int64 = 500)
    n = nx*ns
    c  = zeros(n)                # consumption
    bp = Vector{Float64}(undef, n)      # savings
    counter = 1
    Wout = Vector{Float64}(undef, length(Win))
    qfunction(x::Float64) = mollifier_hetdsge(x, zhi, zlo) #/sumz
    while dist>tol && counter<maxit # for debugging
        # compute c(w) given guess for Win = β*R*E[u'(c_{t+1})]
        for iss in 1:ns
            for ia in 1:nx
                c[nx*(iss-1)+ia] =  min(1/Win[nx*(iss-1)+ia],xgrid[ia]+η)
            end
        end
        bp = repeat(xgrid,ns) - c  # compute bp(w) given guess for Win
        Wout = parameterized_expectations_hetdsge(nx,ns,β,R,ω,H,T,γ,qfunction,xgrid,sgrid,xswts,c,bp,f)
        dist = maximum(abs.(Wout-Win))
        Win = damp*Wout + (1.0-damp)*Win
        counter += 1
    end
    if counter == maxit
        @warn "Euler iteration did not converge"
    end
    tr = kolmogorov_fwd_hetdsge(nx,ns,ω,H,T,R,γ,qfunction,xgrid,sgrid,bp,f)
    return c, bp, Wout, tr
end

function parameterized_expectations_hetdsge(nx::Int,ns::Int, β::AbstractFloat, R::AbstractFloat, ω::AbstractFloat,
                                            H::AbstractFloat, T::AbstractFloat, γ::AbstractFloat,
                                            qfunc::Function,
                                            xgrid::Vector{Float64}, sgrid::Vector{Float64},
                                            xswts::Vector{Float64}, c::Vector{Float64},
                                            bp::Vector{Float64}, f::Array{Float64,2})
    l_out = zeros(nx*ns)
    for iss=1:ns
        for ia=1:nx
            sumn = 0.0
            for isp=1:ns
                for iap=1:nx
                    sumn += (xswts[nx*(isp-1)+iap]/c[nx*(isp-1)+iap])*qfunc((xgrid[iap] - R*(exp(-γ))*bp[nx*(iss-1)+ia] - T)/(ω*H*sgrid[isp]))*f[iss,isp]./sgrid[isp]
                end
            end
            l_out[nx*(iss-1)+ia] = (β*R*(exp(-γ))/ω*H)*sumn
        end
    end
    return l_out
end

function kolmogorov_fwd_hetdsge(nx::Int, ns::Int, ω::AbstractFloat,
                                H::AbstractFloat, T::AbstractFloat,
                                R::AbstractFloat, γ::AbstractFloat,
                                qfunc::Function,
                                xgrid::Vector{Float64},
                                sgrid::Vector{Float64}, bp::Vector{Float64}, f::Array{Float64,2})
    tr = zeros(nx*ns,nx*ns)
    for iss=1:ns
        for ia=1:nx
            for isp=1:ns
                for iap=1:nx
                    tr[nx*(isp-1)+iap,nx*(iss-1)+ia] = qfunc((xgrid[iap] - R*(exp(-γ))*bp[nx*(iss-1)+ia] - T)/(ω*H*sgrid[isp]))*f[iss,isp]./(ω*H*sgrid[isp])
                end
            end
        end
    end
    return tr
end

function mollifier_hetdsge(z::AbstractFloat,ehi::AbstractFloat,elo::AbstractFloat)
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

function dmollifier_hetdsge(x::AbstractFloat, ehi::AbstractFloat, elo::AbstractFloat)
    In = 0.443993816237631
    if x<ehi && x>elo
        temp = (-1.0 + 2.0*(x-elo)/(ehi-elo))
        out  = -(2*temp./((1 - temp.^2).^2)).*(2/(ehi-elo)).*mollifier_hetdsge(x, ehi, elo)
    else
        out = 0.0
    end
    return out
end
