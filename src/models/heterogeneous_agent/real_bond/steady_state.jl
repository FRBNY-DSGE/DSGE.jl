function steadystate!(m::RealBond;
                      βlo::Float64 = 0.4/m[:R],
                      βhi::Float64 = 1 ./m[:R],
                      excess::Float64 = 5000.,
                      tol::Float64 = 1e-5,
                      maxit::Int64 = 100)

    # Load in parameters and grids
    R     = m[:R]
    γ     = m[:γ]
    ν     = m[:ν]
    abar  = m[:abar]
    sigs  = m[:σ_s]

    elo   = get_setting(m, :elo)
    ehi   = get_setting(m, :ehi)

    xgrid = m.grids[:xgrid].points
    xwts  = m.grids[:xgrid].weights
    sgrid = m.grids[:sgrid].points
    swts  = m.grids[:sgrid].weights
    ggrid = m.grids[:ggrid]

    xgrid_total   = m.grids[:xgrid_total]
    sgrid_total   = m.grids[:sgrid_total]
    weights_total = m.grids[:weights_total]

    nx      = get_setting(m, :nx)
    ns      = get_setting(m, :ns)
    n       = get_setting(m, :n)

    # Flag allows one not to recompute steady state
    if get_setting(m, :use_last_βstar)
       βlo = βhi = m[:βstar]
    end

    # Initialization
    count   = 1
    c = zeros(n)
    η = zeros(n)
    ap = zeros(n)
    KF = zeros(n,n)
    μ = zeros(n)
    β = (βlo + βhi)/2.0

    # Define the mollifier
    qfunction(x) = mollifier_realbond(x, ehi, elo)

    # compute constrained consumption and labor supply once and for all
    χss = zeros(nx*ns)
    # ηss = zeros(nx*ns)
    aborrow  = abar/R

    laborzero(s::AbstractFloat, x::AbstractFloat, c::AbstractFloat,
              ν::AbstractFloat, γ::AbstractFloat, aborrow::AbstractFloat) = s^(1+ν) - ((aborrow+c-x)^ν)*c^γ
    chifun(s::AbstractFloat, x::AbstractFloat, ν::AbstractFloat,
           γ::AbstractFloat, aborrow::AbstractFloat) = fzero(c ->laborzero(s, x, c, ν, γ, aborrow), x-aborrow+1e-8, 15.0)
    # etafun(s::AbstractFloat,x::AbstractFloat,ν::AbstractFloat,γ::AbstractFloat,aborrow::AbstractFloat) = chifun(s,x,ν,γ,aborrow)^(-(γ/ν))*s^(1.0/ν)

    for iss in 1:ns
        for ix in 1:nx
            χss[ix + nx*(iss-1)] = chifun(sgrid[iss],xgrid[ix],ν.value,γ.value,aborrow)
            # ηss[ix + nx*(iss-1)] = etafun(sgrid[iss],xgrid[ix],ν,γ,aborrow)
        end
    end

    # guess for the W function W(w) = β R E u'(c_{t+1})
    βguess = 0.8

    # MDχ: This l_in doesn't seem to be used, since it gets immediately overwritten by
    # l_in from policy, unless the while loop condition is never entered
    l_in_const = βguess*R*sum((qfunction.(abar .- xgrid_total).*χss.^(-γ)).*(kron((swts.*ggrid),xwts)))
    l_in = l_in_const*ones(nx*ns)

    while abs(excess) > tol && count < maxit # clearing markets
#        @warn "You're recomputing steady state! BAD!"
        l_in_guess = ones(n)
        β = (βlo + βhi)/2.0

        (c, η, ap, l_in, KF) = policy_realbond(β, R.value, γ.value, ν.value, sigs.value,
                                               ehi, elo, xgrid_total, sgrid_total,
                                               xwts, swts, l_in_guess, ggrid, χss)

        # sspolicy returns the consumption and hours decision rules, a prime
        # decision rule, l_in = updated β R u'(c_{t+1}),
        # KF is the Kolmogorov foward operator
        # and is the map which moves you from cash on hand distribution
        # today to cash on had dist tomorrow

        LPMKF = xwts[1]*swts[1]*KF

        # find eigenvalue closest to 1
        (Dee,Vee) = eigen(LPMKF)
        if abs(Dee[1]-1)>2e-1 # that's the tolerance we are allowing
            @warn "your eigenvalue is ", Dee[1], " which is too far from 1, something is wrong"
        end

        μ = real(Vee[:,1]) # Pick the eigen vector associated with the largest
                           # eigenvalue and move it back to values
                           # mdχ: the μ chosen is the eigenvector associated to the largest
                           # eigenvalue of the KF because that means, we have induced the
                           # shift in μ by the largest amount possible for a given KF.
                           # Hence, to find the steady-state μ, we want to pick the
                           # potential μ that changes by the most (i.e. the largest
                           # eigenvector) so we have an upper bound on the degree of change
                           # of μ induced by the KF.

        μ = μ/(weights_total'*μ) # Scale of eigenvectors not determinate: rescale to integrate to exactly 1

        excess = weights_total'*(μ.*ap)  # compute excess supply of savings, which is a fn of w

        # bisection
        if excess > 0
            βhi = β
        elseif excess < 0
            βlo = β
        end
        count += 1
    end

    m[:lstar] = l_in
    m[:cstar] = c
    m[:ηstar] = η
    m[:μstar] = μ
    m[:βstar] = β
    m[:χstar] = χss

    nothing
end

function policy_realbond(β::AbstractFloat, R::AbstractFloat,
                         γ::AbstractFloat, ν::AbstractFloat,
                         sigs::AbstractFloat, ehi::AbstractFloat, elo::AbstractFloat,
                         xgrid_total::Vector{Float64}, sgrid_total::Vector{Float64},
                         xwts::Vector{Float64}, swts::Vector{Float64},
                         l_in::Vector{Float64}, ggrid::Vector{Float64},
                         χss::Vector{Float64};
                         damp::Float64 = 0.4, dist::Float64 = 1.,
                         tol::Float64 = 1e-4, maxit::Int64 = 500)
    # outputs: [c, η, ap, l_out, tr]

    nx      = length(xwts)
    ns      = length(swts)
    n       = ns*nx
    c       = zeros(n)      # consumption
    η       = zeros(n)      # hours
    ap      = zeros(n)      # savings
    counter = 1
    l_in    = ones(n)
    l_out   = copy(l_in)    # initialize l_out
    qfunction(x) = mollifier_realbond(x,ehi,elo)

    while dist > tol && counter < maxit
        # compute c(w) given guess for l_in = β*R*E[u'(c_{t+1})]

        for i in 1:n
            c[i]  = min(l_in[i]^(-1.0/γ), χss[i])
        end

        η  = (sgrid_total.^(1.0/ν)).*(c.^(-γ/ν))
        ap = R*(xgrid_total+sgrid_total.*η - c)

        l_out = parameterized_expectations_realbond(nx, ns, qfunction, xgrid_total, ggrid, xwts, swts, ap, c, γ, β, R)

        dist = norm(l_out - l_in, Inf)

        # if counter%1==0
        #     println(dist)
        # end
        l_in = damp*l_out + (1.0-damp)*l_in
        counter += 1
    end

    if counter == maxit
        @warn "Euler iteration did not converge"
    end

    # c  = min.(l_in.^(-1.0/γ),χss)
    # η  = (sgrid_total.^(1.0/ν)).*(c.^(-γ/ν))
    # ap = R*(xgrid_total+sgrid_total.*η - c)

    tr = kolmogorov_fwd_policy_realbond(nx, ns, qfunction, xgrid_total, ggrid, ap)

    return c, η, ap, l_out, tr
end

function parameterized_expectations_realbond(nx::Int, ns::Int,
                                             qfunction::Function,
                                             xgrid_total::Vector{Float64}, ggrid::Vector{Float64},
                                             xwts::Vector{Float64}, swts::Vector{Float64},
                                             ap::Vector{Float64}, c::Vector{Float64},
                                             γ::Float64, β::Float64, R::Float64)
    l_out   = ones(ns*nx)    # initialize l_out
    for ix in 1:nx
        for iss in 1:ns
            sumn = 0.0
            for ixp in 1:nx
                for isp in 1:ns
                    # This is the parameterized expectations algorithm approximation of the
                    # conditional expectation within the Euler equation
                   sumn += qfunction(ap[ix+nx*(iss-1)] - xgrid_total[nx*(isp-1)+ixp])*ggrid[isp]*xwts[ixp]*swts[isp]*(c[nx*(isp-1)+ixp]^(-γ))
                end
            end
            l_out[ix+nx*(iss-1)] = β*R*sumn
        end
    end

    return l_out
end

function kolmogorov_fwd_policy_realbond(nx::Int, ns::Int,
                                        qfunction::Function,
                                        xgrid_total::Vector{Float64}, ggrid::Vector{Float64},
                                        ap::Vector{Float64})
    tr = zeros(ns*nx,ns*nx)    # transition matrix (KFE) mapping todays dist of w to w' (i.e. μ to μ′)
    for ix in 1:nx
        for ixp in 1:nx
            for iss in 1:ns
                for isp in 1:ns
                    tr[ixp+nx*(isp-1),ix+nx*(iss-1)] += qfunction(ap[ix+nx*(iss-1)] - xgrid_total[nx*(isp-1)+ixp])*ggrid[isp]
                end
            end
        end
    end
    return tr
end

function mollifier_realbond(z::AbstractFloat,ehi::AbstractFloat,elo::AbstractFloat)
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

function dmollifier_realbond(x::AbstractFloat, ehi::AbstractFloat, elo::AbstractFloat)
    In = 0.443993816237631
    if x<ehi && x>elo
        temp = (-1.0 + 2.0*(x-elo)/(ehi-elo))
        out  = -(2*temp./((1 - temp.^2).^2)).*(2/(ehi-elo)).*mollifier_realbond(x, ehi, elo)
    else
        out = 0.0
    end
    return out
end
