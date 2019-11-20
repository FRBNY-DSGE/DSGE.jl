@inline function policy(m::KrusellSmith,
                        net_return::Float64, wage::Float64,
                        l_in::Vector{Float64}, g::Vector{Float64}, KF::Array{Float64,2};
                        damp::Float64 = 0.8, tol::Float64 = 1e-5)
    @assert 0 <= damp <= 1

    nw::Int64                = get_setting(m, :nw)
    ns::Int64                = get_setting(m, :ns)
    wgrid::DSGE.Grid         = get_grid(m, :wgrid)
    sgrid::DSGE.Grid         = get_grid(m, :sgrid)

    zhi::Float64             = get_setting(m, :zhi)
    zlo::Float64             = get_setting(m, :zlo)
    smoother::Float64        = get_setting(m, :smoother)
    γ::Float64               = m[:γ]
    β::Float64               = m[:β]

    dist::Float64            = 1.0         # distance between l_in and l_out on decision rules
    l_out::Array{Float64,1}  = copy(l_in)  # initialize l_out

    ap = Array{Float64,1}(undef, nw) # savings
    q  = Array{Float64,2}(undef, nw, nw * ns) # auxilary variable to compute LHS of euler

    count = 1

    while dist > tol && count < 5000

        # compute c(w) given guess for l_in = β*R*E[u'(c_{t+1})]
        c::Array{Float64,1} = min.(l_in.^(-1 / γ), wgrid.points)

        # compute ap(w) given guess for l_in
        ap = wgrid.points - c

        # The term contained inside the (⋅) next to g in equation (9)
        q = euler_lhs(q, ap, wgrid.points, sgrid.points, zhi, zlo, smoother, wage, net_return, nw, ns)

        # Equation (9)
        l_out = β * net_return * q * kron(sgrid.weights .* (g / wage),
                                              wgrid.weights .* (c.^-γ))
        dist = maximum(abs.(l_out - l_in))
        l_in = damp * l_out + (1.0 - damp) * l_in;

        count += 1
    end
    if count == 5000
        @warn "Euler iteration did not converge"
    end

    # Equation (10)
    # You can re-use most of the functional form of the Euler equation to
    # back out the KF for μ(w)
    KF = kolmogorov_forward(q, KF, g, sgrid.weights, wage, nw, ns)

    return ap, l_out, KF
end

function steadystate!(m::KrusellSmith, tol::Float64 = 1e-5, maxit::Int64 = 100, Kdamp::Float64 = 0.1)
    # Load settings and grids
    nw::Int64  = get_setting(m, :nw)

    wgrid::DSGE.Grid = get_grid(m, :wgrid)
    sgrid::DSGE.Grid = get_grid(m, :sgrid)

    zhi::Float64      = get_setting(m, :zhi)
    zlo::Float64      = get_setting(m, :zlo)
    smoother::Float64 = get_setting(m, :smoother)

    α::Float64 = m[:α]
    δ::Float64 = m[:δ]
    β::Float64 = m[:β]

    g = trunc_lognpdf(sgrid, m[:μ_s].value, m[:σ_s].value)      # Density g (of skill, s) evaluated on sgrid

    m[:Lstar] = L::Float64 = quadrature_sum(g, sgrid)           # Aggregate Labor endowment

    K::Float64 = L*(α /((1.0 / β) - 1.0 + δ))^(1.0 / (1.0 - α)) # Aggregate Capital Stock

    l_in = fill(0.3, nw)        # This corresponds to l^i_t = βE_t[R_{t+1} u'(c^i_{t+1})]
                                #                           = βE_t[R_{t+1} min(l^i_{t+1}, (w^i_{t+1}^(-1/γ))]

    ap   = Array{Float64,1}(undef, nw)
    KF   = Array{Float64,2}(undef, nw, nw) # The Kolmogorov Forward Equation for μ
    μ    = Array{Float64,1}(undef, nw)     # The cross-sectional density of cash on hand, w
    newK = 0.0

    diseqm = 1000.
    count = 1

    # Equation (11)
    while abs(diseqm) > tol && count < maxit # clearing markets

        wage = MPL(1.0, L, K, α)
        net_return = MPK(1.0, L, K, α, δ)

        ap, l_in, KF = policy(m, net_return, wage, l_in, g, KF)

        # sspolicy returns the consumption decision rule, a prime
        # decision rule, l_in = updated β R u'(c_{t+1}),
        # KF is the Kolmogorov forward operator and is the map which moves you
        # from cash on hand distribution today to cash on hand dist tomorrow

        # find eigenvalue closest to 1
        D,V = eigen(wgrid.weights[1] * KF)

        if norm(D[1] - 1.0) > 2e-1 # that's the tolerance we are allowing
            @warn "your eigenvalue is too far from 1, something is wrong"
        end
        μ = real(@view(V[:,1]))  # Pick the eigenvector associated with the largest
                                 # eigenvalue and moving it back to values

        μ = μ / (wgrid.weights' * μ)        # Scale of eigenvectors not determinate:
                                            # rescale to integrate to exactly 1
        newK  = wgrid.weights' * (μ .* ap)  # compute excess supply of savings, which is a fn of w
        diseqm = newK - K

        K += Kdamp * diseqm
        count += 1
    end

    m[:lstar] = l_in
    m[:cstar] = min.(m[:lstar].value.^(-1.0/m[:γ]), wgrid.points)
    m[:μstar] = μ
    m[:Kstar] = newK
    m[:KFstar] = KF

    nothing
end

@inline function mollifier(z::Float64,zhi::Float64,
                   zlo::Float64, smoother::Float64; In::Float64=0.443993816237631)
    # mollifier function
    if z < zhi && z > zlo
        temp = (-1.0 + 2.0 * (z - zlo) / (zhi - zlo)) / smoother
        out = (2.0/(zhi - zlo)) * exp(-1.0 / (1.0 - temp*temp)) / (In * smoother);
    else
        out = 0.0
    end
    return out
end

# Derivative of the mollifier
@inline function dmollifier(x::Float64, zhi::Float64, zlo::Float64, smoother::Float64;
                            In::Float64=0.443993816237631)
    temp = (-1.0 + 2.0*(x-zlo)/(zhi-zlo))/smoother
    dy  = -(2*temp./((1 - temp.^2).^2)).*(2/(smoother*(zhi-zlo))).*mollifier(x, zhi, zlo, smoother)
end

@inline function trunc_lognpdf(grid::Grid, mu::Float64, sig::Float64)
    x = grid.points
    xlo = x[1]
    xhi = x[end]

    # truncated log normal pdf
    logn = LogNormal(mu, sig)
    return (x .- xlo .>= -eps()) .* (x .- xhi .<= eps()) .* pdf.(logn,x) / (cdf(logn,xhi) .- cdf(logn,xlo))
end

# ω_t
@inline function MPL(Z::Float64, L::Float64, K::Float64, α::Float64)
    (1.0 - α) * Z * (K^α) * (L^(-α))  # Marginal product of labor
end

# r_t + (1 - δ)
@inline function MPK(Z::Float64, L::Float64, K::Float64, α::Float64, δ::Float64)
    α * Z * (K^(α - 1.0)) * (L^(1.0 - α)) + (1.0 - δ)   # Marginal product of capital + 1-δ
end

####### Is this function ever used? #########
#=
# Re-computing
@inline function kolmogorov_forward(g::Vector{Float64}, sweights::Vector{Float64},
                            ap::Vector{Float64}, wgrid::Vector{Float64}, sgrid::Vector{Float64},
                            zhi::Float64, zlo::Float64, smoother::Float64,
                            wage::Float64, R::Float64, nw::Int64, ns::Int64)

    KF = Array{Float64, 2}(nw, nw)   # transition matrix mapping todays dist of w to w'
    for iw=1:nw
        for iwp=1:nw
            sumns = 0.0
            cursum = (wgrid[iwp] - R*ap[iw]) / wage
            for isp=1:ns
                sumns += mollifier(cursum - sgrid[isp], zhi, zlo, smoother) * (g[isp]) * sweights[isp]
            end
            @inbounds KF[iwp,iw] = sumns / wage
        end
    end
    return KF
end
=#

# Indexing into q since its already been computed
@inline function kolmogorov_forward(q::Matrix{Float64}, KF::Array{Float64,2}, g::Vector{Float64},
                                    sweights::Vector{Float64}, wage::Float64, nw::Int64, ns::Int64)
    g_sweights_prod::Array{Float64,1} = g / wage .* sweights
    nw_array = nw * (0:ns-1)
    for iw = 1:nw
        for iwp = 1:nw
            @inbounds KF[iwp, iw] = sum(q[iw, nw_array .+ iwp] .* g_sweights_prod)
        end
    end
    return KF
end

@inline function euler_lhs(q::Matrix{Float64}, ap::Vector{Float64}, wgrid::Vector{Float64},
                   sgrid::Vector{Float64}, zhi::Float64, zlo::Float64, smoother::Float64,
                   wage::Float64, net_return::Float64, nw::Int64, ns::Int64)
    for iwp = 1:nw
        for iw = 1:nw
            wgrid_term = (wgrid[iwp] - net_return * ap[iw]) / wage
            for iss = 1:ns
                @inbounds q[iw, nw * (iss - 1) + iwp] = mollifier(wgrid_term - sgrid[iss], zhi, zlo, smoother)
            end
        end
    end
    return q
end
