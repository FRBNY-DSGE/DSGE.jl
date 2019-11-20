function steadystate!(m::HetDSGE;
                      βlo::Float64 = 0.5*exp(m[:γ])/(1 + m[:r]),
                      βhi::Float64 = exp(m[:γ])/(1 + m[:r]),
                      excess::Float64 = 5000.,
                      tol::Float64 = 1e-4,
                      maxit::Int64 = 20,
                      βband::Float64 = 1e-2)

    # Load settings
    nx = get_setting(m, :nx)
    ns = get_setting(m, :ns)

    xlo = get_setting(m, :xlo)
    xhi = get_setting(m, :xhi)
    zlo = get_setting(m, :zlo)
    zhi = get_setting(m, :zhi)

    # Load parameters
    R = 1 + m[:r].value
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

    βlo_temp = βlo
    βhi_temp = βhi

    counter = 1
    n  = ns*nx
    c  = zeros(n)
    bp = zeros(n)
    KF = zeros(n,n)
    μ  = zeros(n)

    # Initial guess
    β   = 1.0
    Win = Vector{Float64}(undef, n)
    Win_guess = ones(n)

    if get_setting(m, :use_last_βstar) && !isnan(m[:βstar].value)
        βlo = βhi = m[:βstar].value
    elseif !isnan(m[:βstar].value)
        βlo_temp = m[:βstar].value - βband
        βhi_temp = m[:βstar].value + βband

        c, bp, Win, KF = policy_hetdsge(nx, ns, βlo_temp, R, ω, H, η, T, γ,
                                        zhi, zlo, xgrid, sgrid, xswts, Win_guess, f)
        excess_lo, μ = compute_excess(xswts, KF, bp)

        if excess_lo < 0 && abs(excess_lo) > tol
            βlo = βlo_temp
            c, bp, Win, KF = policy_hetdsge(nx, ns, βhi_temp, R, ω, H, η, T, γ, zhi,
                                            zlo, xgrid, sgrid, xswts, Win_guess, f)
            excess_hi, μ = compute_excess(xswts, KF, bp)

            if excess_hi > 0
                βhi = βhi_temp
            end
        elseif excess_lo >= 0
            βhi = βlo_temp
        end
    end

    while abs(excess) > tol && counter < maxit # clearing markets
        β = (βlo + βhi) / 2.0
        c, bp, Win, KF = policy_hetdsge(nx, ns, β, R, ω, H, η, T, γ, zhi, zlo, xgrid, sgrid,
                                        xswts, Win_guess, f)
        excess, μ = compute_excess(xswts, KF, bp)

        # bisection
        if excess > 0
            βhi = β
        elseif excess < 0
            βlo = β
        end
        counter += 1
        if get_setting(m, :use_last_βstar) && !isnan(m[:βstar].value)
            break
        end
    end
    m[:lstar] = Win
    m[:cstar] = c
    m[:μstar] = μ
    m[:βstar] = β
    m[:β_save] = β

    #return Win, c, μ, β
    nothing
end

@inline function compute_excess(xswts::Vector{T}, KF::Matrix{T}, bp::Vector{T}) where T<:Float64
    LPMKF = xswts[1] * KF
    # Find eigenvalue closest to 1
    (D,V) = (eigen(LPMKF)...,)
    max_D = argmax(abs.(D))
    D = D[max_D]

    if abs(D - 1) > 2e-1 # that's the tolerance we are allowing
        @warn "Your eigenvalue is too far from 1, something is wrong."
    end
    # Pick eigenvector associated w/ largest eigenvalue and moving it back to values
    μ = real(V[:,max_D])
    μ = μ ./ dot(xswts, μ) # Scale of eigenvectors not determinate: rescale to integrate to 1
    excess = dot(xswts, (μ .* bp)) # compute excess supply of savings, which is a fn of w

    return excess, μ
end

function policy_hetdsge(nx::Int, ns::Int, β::AbstractFloat, R::AbstractFloat,
                        ω::AbstractFloat, H::AbstractFloat,
                        η::AbstractFloat, T::AbstractFloat,
                        γ::AbstractFloat, zhi::AbstractFloat,
                        zlo::AbstractFloat,
                        xgrid::Vector{Float64}, sgrid::Vector{Float64},
                        xswts::Vector{Float64}, Win::Vector{Float64},
                        f::Array{Float64,2}, damp::Float64 = 0.5, dist::Float64 = 1.,
                        tol::Float64 = 1e-4, maxit::Int64 = 500)
    n    = nx*ns
    c    = zeros(n)                # consumption
    bp   = Vector{Float64}(undef, n)      # savings
    Wout = Vector{Float64}(undef, length(Win))
    counter = 1
    qfunction_hetdsge(x::Float64) = mollifier_hetdsge(x, zhi, zlo) #/sumz

    while dist>tol && counter<maxit # for debugging
        # compute c(w) given guess for Win = β*R*E[u'(c_{t+1})]
        for iss in 1:ns
            for ia in 1:nx
                c[nx*(iss-1)+ia] = min(1/Win[nx*(iss-1)+ia],xgrid[ia]+η)
            end
        end
        bp = repeat(xgrid, ns) - c  # compute bp(w) given guess for Win
        Wout = parameterized_expectations_hetdsge(nx, ns, β, R, ω, H, T, γ,
                                                  qfunction_hetdsge, xgrid, sgrid, xswts, c, bp, f)

        dist = maximum(abs.(Wout - Win))
        Win  = damp*Wout + (1.0-damp) * Win
        counter += 1
    end
    if counter == maxit
        @warn "Euler iteration did not converge"
    end
    tr = kolmogorov_fwd_hetdsge(nx, ns, ω, H, T, R, γ, qfunction_hetdsge, xgrid, sgrid, bp, f)
    return c, bp, Wout, tr
end

@inline function parameterized_expectations_hetdsge(nx::Int, ns::Int, β::AbstractFloat,
                                                    R::AbstractFloat, ω::AbstractFloat,
                                                    H::AbstractFloat, T::AbstractFloat,
                                                    γ::AbstractFloat,
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
                    sumn += (xswts[nx*(isp-1)+iap]/c[nx*(isp-1)+iap]) *
                        qfunc((xgrid[iap] - R*(exp(-γ))*bp[nx*(iss-1)+ia] - T) /
                              (ω*H*sgrid[isp])) * f[iss,isp] ./ sgrid[isp]
                end
            end
            l_out[nx*(iss-1)+ia] = (β*R*(exp(-γ))/ω*H)*sumn
        end
    end
    return l_out
end

@inline function kolmogorov_fwd_hetdsge(nx::Int, ns::Int, ω::AbstractFloat,
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
                    tr[nx*(isp-1)+iap, nx*(iss-1)+ia] = qfunc((xgrid[iap] - R*(exp(-γ))*bp[nx*(iss-1)+ia] - T)/(ω*H*sgrid[isp])) * f[iss,isp] ./ (ω * H * sgrid[isp])
                end
            end
        end
    end
    return tr
end

@inline function mollifier_hetdsge(z::AbstractFloat,ehi::AbstractFloat,elo::AbstractFloat)
    # mollifier function
    In = 0.443993816237631
    if z<ehi && z>elo
        temp = -1.0 + 2.0 * (z - elo) / (ehi - elo)
        return (2.0 / (ehi - elo)) * exp(-1.0 / (1.0 - temp^2)) / In
    end
    return 0.0
end

@inline function dmollifier_hetdsge(x::AbstractFloat, ehi::AbstractFloat, elo::AbstractFloat)
    In = 0.443993816237631
    if x<ehi && x>elo
        temp = (-1.0 + 2.0*(x-elo)/(ehi-elo))
        out  = -(2*temp./((1 - temp.^2).^2)).*(2/(ehi-elo)).*mollifier_hetdsge(x, ehi, elo)
    else
        out = 0.0
    end
    return out
end
