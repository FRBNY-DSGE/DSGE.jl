"""
Sparse identity matrix - since deprecated in 0.7
"""
function my_speye(n::Integer)
    return SparseMatrixCSC{Float64}(LinearAlgebra.I, n, n)
end

"""
Sparse identity matrix - since deprecated in 0.7
"""
function my_speye(T::Type, n::Integer)
    return SparseMatrixCSC{T}(LinearAlgebra.I, n, n)
end

function stat_dist(λ::AbstractArray{T}) where {T<:Float64}
    F     = eigen(copy(λ))
    eigvl = F.values
    eigvc = F.vectors
    pos   = argmin(abs.(eigvl))
    statdist = eigvc[:,pos]
    return statdist ./ sum(statdist)
end

function changeBasis(basis,inv_basis,g1,psi,pi,c,g0)
    # Update dynamics with new basis
    # by SeHyoun Ahn, June 2016

    g1 = basis*g1*inv_basis
    pi = basis*pi
    psi = basis*psi
    c = basis*c
    if (nargin==7)
        g0 = basis*g0*inv_basis
    else
        g0 = my_speye(size(g1, 1))
    end
    return g1, psi, pi, c, g0
end

function adj_cost_fn(d, a_grid, χ0, χ1, χ2, a_lb)
    d_scaled = abs.(d ./ max.(a_grid, a_lb))
    adj_cost = max.(a_grid, a_lb) .* (χ0 * (d_scaled) .+ 1.0 / (1.0 .+ χ2) *
                                    ((d_scaled).^(1. +χ2) * χ1.^(-χ2)))
    return adj_cost
end

function cleanG0Sparse(g0,g1,c,pi,psi)
    # Inverts g0. For certain form on non-invertible g0 coming from constraints,
    # it will compute the proper inverse.
    #
    # by SeHyoun, June 2016

    n = size(g0,1)

    tmp = (max(abs(g0),[],2)==0)
    #tmp=(max(abs([g0,psi]),[],2)==0)

    redundant = findall(tmp)
    keep=findall(1 .- tmp)
    n_keep = length(keep)

    base = sparse(n,n_keep)
    base[keep,:] = my_speye(n_keep)

    base[redundant,:] = -g1[redundant,redundant] \ g1[redundant,keep]

    inv_base = sparse(n_keep,n)
    inv_base[:,keep] = my_speye(n_keep)

    g0=inv_base*g0*base
    g1=inv_base*g1*base
    g1=g0\g1
    psi=g0\inv_base*psi
    pi=g0\inv_base*pi
    c=g0\inv_base*c
    return base, inv_base, g1, c, pi, psi
end

function  covv(x,y)
    cov_xy = cov(x,y)
    cov_xy = cov_xy[1,2]
    return cov_xy
end

function Gini_fn(gg,c,dab_aux)

    S_c  = cumsum(gg .* dab_aux .* c[:])/sum(gg .* dab_aux .* c[:])
    trapez_c = 0.5*(S_c[1]*gg[1]*dab_aux[1] +
                      sum((S_c[2:end] + S_c[1:end-1]).*gg[2:end].*dab_aux[2:end]))
    C_Gini = 1 - 2*trapez_c
    return C_Gini
end

function gstate(G1, impact; pick = Matrix{Float64}(undef, 0, 0))
    #function  [GS,pickn,ok,uis,vs]=gstate(G1,impact,pick)
    # G1:    the coefficient on  lagged y from the output of gensys.m.
    # impact:the coefficient on exogenous shocks from the output of gensys.m
    # pick:  an optional guess at a matrix of coefficients that extracts
    #        the state vector from y
    # ok:     0 if pick matrix is just no good.
    #         1 if pick matrix is usable forward, but loses initial date information
    #         2 if pick matrix is not usable forward, is part of a correct state vector,
    #           but is not complete
    #         3 if pick matrix is     usable forward, is part of a correct state vector,
    #           but is not complete
    #         4 if pick matrix is not usable forward, summarizes past, but is redundant
    #         5 if pick matrix is     usable forward, summarizes past, but is redundant
    #         6 if pick matrix summarizes past and is not redundant, but is not usable forward
    #         7 if pick matrix is perfect, both forward and as history summary
    # pickn: a matrix of coefficients that extracts the state vector from y.  Equal
    #        to pick if pick is supplied and ok=1; otherwise, an ok-maximizing pick matrix that
    #        is somewhat similar to pick.
    # GS:    the matrix of coefficients on the lagged state variable
    # uis,vs: If uis'*vs is square and full rank, ok=7 is possible, otherwise not.
    #         If ok<2, an ok>2 can be found by trying a pick in the row space of vs'.
    #         Any pick with pick*uis full column rank will provide a foward state
    #         (i.e. ok an odd number).
    # uiu:   uiu'y(t)=0, and this is interpretable as the decision rule setting controls
    #        as functions of the states.
    # The solution was in the form y(t)=G1*y(t-1)+impact*z(t).
    # Now it's in the form y(t)=GS*pickn*y(t-1)+impact*z(t).
    # In general, pickn is mxn with m<n.
    REALSMALL = 1e-9
    nr, nc = size(G1)
    F = svd(G1)
    u, d, v = F.U, F.S, F.V #[u d v] = svd(G1)
    top = findall(d .> REALSMALL) #diag(d)>REALSMALL)
    nd = top[end]
    us = u[:,top]
    uu = u[:,top+1:end]
    d  = d[top,top]
    vs = v[:,top]
    if length(pick) == 0 #nargin <= 2
        pick = vs'
        vp   = v
        vps  = vs
        dp   = Matrix{Float64}(I,nd)
        ups  = Matrix{Float64}(I,nr,nd)
        ndp  = nd
    else
        Fp = svd(pick)
        up, dp, vp = Fp.U, Fp.S, Fp.V #[up dp vp] = svd(pick)
        topp = findall(dp .> REALSMALL)
        ndp  = topp[end]
        vps  = vp[:,topp]
        dp   = dp[topp,topp]
        ups  = up[:,topp]
    end

    # If we were worried about efficiency, we'd skip some of this when pick = vs.
    #Does pick summarize history?  (pick in v' row space, v in vp' row space).
    pinv = all(all((pick'-vs*vs'*pick').^2<REALSMALL))
    vinp = all(all((vs-vps*vps'*vs).^2<REALSMALL))
    okpast = pinv+vinp
    #  Does pick summarize all current info?  (impact in us column space, pick*uu full rank)
    Fi = svd(hcat(impact, us)) # RECA: [impact us]
    ui, di, vi = Fi.U, Fi.S, Fi.V
    topi = findall(di .> REALSMALL)
    ndi = length(topi)
    uis = ui[:,topi]
    uiu = ui[:,ndi+1:size(ui,2)]
    if ndi < size(G1,1)
        if size(pick,1) < size(uis,2)
            oknow = 0
        else
            Ft = svd(pick*uis)
            ut, dt, vt = Ft.U, Ft.S, Ft.V
            toppu = findall(dt .> REALSMALL)
            if length(toppu) < size(us,2)
                oknow = 0
            else
                oknow = 1
            end
        end
    else
        oknow = 0
    end
    if vinp
        GS = G1/pick
        pickn = pick
    elseif pinv
        r = vs-vps*vps'*vs
        Fr = svd(r)
        ur, dr, vr = Fr.U, Fr.S, Fr.V
        topr = findall(dr .> REALSMALL)
        p2 = ur[:,topr]
        pickn = [pick;p2']
        GS = G1/pickn
    elseif oknow
        GS = G1/[pick;uiu']
        GS = GS[:,1:size(pick,1)]
        pickn = pick
    else
        pickn = vs'
        GS = us*d
    end
    ok = oknow+2*pinv+4*vinp
    return GS,pickn,ok,uis,uiu,vs
end

function opt_deposits(Va, Vb, a_grid, χ0, χ1, χ2, a_lb)
    indx_plus  = ((Va ./ Vb .- 1 .- χ0) .> 0)
    indx_minus = ((Va ./ Vb .- 1 .+ χ0) .< 0)

    return χ1 .* (max.(Va ./ Vb .- 1 .- χ0, 0)) .^ (1/χ2) .*
        max.(a_grid, a_lb) .* indx_plus .+ (-χ1) .*
        (max.(-(Va./Vb .- 1) .- χ0, 0)) .^ (1/χ2) .* max.(a_grid, a_lb) .* indx_minus
end


function opt_deposits(Va::T, Vb::T, a::T, χ0::T, χ1::T, χ2::T) where {T<:Real}
    indx_plus  = ((Va ./ Vb .- 1 .- χ0) .> 0)
    indx_minus = ((Va ./ Vb .- 1 .+ χ0) .< 0)

    return χ1 .* (max.(Va ./ Vb .- 1 .- χ0, 0)) .^ (1/χ2) .* a .* indx_plus .+
        (-χ1) .* (max.(-(Va./Vb .- 1) .- χ0, 0)) .^ (1/χ2) .* a .* indx_minus
end


@inline function u_fn(c::Array{T}, γ::Float64) where T
    if γ == 1.0
        u = log.(c)
    else
        u = 1.0 / (1.0 - γ) * (c .^ (1-γ) - 1.0)
    end
    return u
end
#=
"""
    <(a::Complex, b::Complex)

Compare real values of complex numbers.
"""
function <(a::Complex, b::Complex)
    return a.re < b.re
end

"""
    <(a::Real, b::Complex)

Compare real values of complex numbers.
"""
function <(a::Real, b::Complex)
    return a < b.re
end

"""
    <(a::Complex, b::Real)

Compare real values of complex numbers.
"""
function <(a::Complex, b::Real)
    return a.re < b
end

function min(a::Complex, b::Real)
    return min(a.re, b)
end

function min(a::Complex, b::Complex)
    return min(a.re, b.re)
end

function min(a::Real, b::Complex)
    return min(a, b.re)
end

function max(a::Complex, b::Real)
    return max(a.re, b)
end

function max(a::Complex, b::Complex)
    return max(a.re, b.re)
end

function max(a::Real, b::Complex)
    return max(a, b.re)
end
=#
