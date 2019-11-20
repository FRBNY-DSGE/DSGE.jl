"""
```
function my_speye(n::Integer)
```
Sparse identity matrix - since deprecated in 0.7
"""
@inline function my_speye(n::Integer)
    return SparseMatrixCSC{Float64}(LinearAlgebra.I, n, n)
end

"""
```
function my_speye(T::Type, n::Integer)
```
Sparse identity matrix - since deprecated in 0.7
"""
@inline function my_speye(T::Type, n::Integer)
    return SparseMatrixCSC{T}(LinearAlgebra.I, n, n)
end

"""
```
function stat_dist(λ::AbstractArray{T}) where {T<:Float64}
```
Computes stationary distribution.
"""
function stat_dist(λ::AbstractArray{T}) where {T<:Float64}
    F     = eigen(copy(λ))
    eigvl = F.values
    eigvc = F.vectors
    pos   = argmin(abs.(eigvl))
    statdist = eigvc[:,pos]
    return statdist ./ sum(statdist)
end

"""
```
@inline function changeBasis(basis,inv_basis,g1,psi,pi,c,g0)
```
Update dynamics with new basis.
(by SeHyoun Ahn, June 2016)

"""
@inline function changeBasis(basis,inv_basis,g1,psi,pi,c,g0)
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

"""
```
function stateSpaceReduction(g0, g1, n_v, n_g, n_p, n_Z, reduceDist_hor)
```
Reduce dimensionality of state space.

# Arguments
- `g0`: LHS matrix, only used to check it satisfies require form
- `g1`: Dynamics matrix
- `n_v::Int64`: number of jump variables
- `n_g::Int64`: number of state variables
- `n_p::Int64`: number of static constraints
- `n_Z::Int64`: number of exogenous variables

# Output
- `state_red`: tranformation to get from full grid to reduced states
- `inv_state_red`: inverse transform
- `n_g`: number of state variables after reduction
"""
function stateSpaceReduction(g0, g1, n_v, n_g, n_p, n_Z, reduceDist_hor)
    # Check to make sure that LHS satisfies the required form
    loc = [1:n_v+n_g, n_v+n_g+n_p+1:n_v+n_g+n_p+n_Z]
    #if (maximum(abs.(g0[loc,loc] .- my_speye(n_v+n_g+n_Z)))≈0)
    #    error("Make sure that g0 is normalized.")
    #end

    n_total = n_v + n_g

    # Slice Dynamics Equation into Different Parts
    B_pv = g1[n_v+n_g+1:n_v+n_g+n_p,n_v+n_g+1:n_v+n_g+n_p] \ g1[n_v+n_g+1:n_v+n_g+n_p,1:n_v]
    B_pZ = g1[n_v+n_g+1:n_v+n_g+n_p,n_v+n_g+1:n_v+n_g+n_p] \ g1[n_v+n_g+1:n_v+n_g+n_p,n_v+n_g+n_p+1:n_v+n_g+n_p+n_Z]
    B_pg = g1[n_v+n_g+1:n_v+n_g+n_p,n_v+n_g+1:n_v+n_g+n_p] \ g1[n_v+n_g+1:n_v+n_g+n_p,n_v+1:n_v+n_g]
    B_gg = g1[n_v+1:n_v+n_g,n_v+1:n_v+n_g]
    B_gp = g1[n_v+1:n_v+n_g,n_v+n_g+1:n_v+n_g+n_p]

    # Orthonormalize B_pg as the direction to keep
    # Drop redundant Directions
    ~, d0, V_g = svd(B_pg)#,'econ')
    aux = diagm(0 => d0)
    n_Bpg = Int64(sum(aux .> 10 * eps() * aux[1]))
    V_g = V_g[:, 1:n_Bpg] .* aux[1:n_Bpg]'

    ## Arnoldi Iteration
    hor = reduceDist_hor
    A(x::Matrix{Float64}) = B_gg' * x - B_pg' * (B_gp'*x)
    # [V_g,~] = block_arnoldi_func(A,V_g,hor)
    V_g, ~ = deflated_block_arnoldi(A, V_g, hor)
    n_g = size(V_g, 2)

    ## Build State-Space Reduction transform
    state_red = spzeros(Float64, n_v + n_g + n_Z, n_total + n_p + n_Z)
    state_red[1:n_v,1:n_v] = my_speye(n_v)
    state_red[n_v+1:n_v+n_g,n_v+1:n_total] = V_g'
    state_red[:,n_total+1:n_total+n_p] = 0
    state_red[n_v+n_g+1:n_v+n_g+n_Z,n_total+n_p+1:n_total+n_p+n_Z] = my_speye(n_Z)

    ## Build inverse transform
    inv_state_red = spzeros(Float64, n_total+n_p+n_Z,n_v+n_g+n_Z)
    inv_state_red[1:n_v,1:n_v] = my_speye(n_v)
    inv_state_red[n_v+1:n_total,n_v+1:n_g+n_v] = V_g
    inv_state_red[n_total+1:n_total+n_p,n_v+1:n_v+n_g] = -B_pg*V_g
    inv_state_red[n_total+1:n_total+n_p,1:n_v] = -B_pv
    inv_state_red[n_total+1:n_total+n_p,n_v+n_g+1:n_v+n_g+n_Z] = -B_pZ
    inv_state_red[n_total+n_p+1:n_total+n_p+n_Z,n_v+n_g+1:n_v+n_g+n_Z] = my_speye(n_Z)

    return state_red, inv_state_red, n_g
end

"""
```
function cleanG0Sparse(g0, g1, c, pi, psi)
```
Inverts g0. For certain form on non-invertible g0 coming from constraints,
it will compute the proper inverse.

(by SeHyoun, June 2016)
"""
function cleanG0Sparse(g0, g1, c, pi, psi)
    n = size(g0,1)

    tmp = (max(abs(g0),[],2)==0)

    redundant = findall(tmp)
    keep=findall(1 .- tmp)
    n_keep = length(keep)

    base = sparse(n,n_keep)
    base[keep,:] = my_speye(n_keep)

    base[redundant,:] = -g1[redundant,redundant] \ g1[redundant,keep]

    inv_base = sparse(n_keep,n)
    inv_base[:,keep] = my_speye(n_keep)

    g0  = inv_base * g0 * base
    g1  = inv_base * g1 * base
    g1  = g0 \ g1
    psi = g0 \ inv_base * psi
    pi  = g0 \ inv_base * pi
    c   = g0 \ inv_base * c

    return base, inv_base, g1, c, pi, psi
end

function covv(x,y)
    cov_xy = cov(x,y)
    cov_xy = cov_xy[1,2]
    return cov_xy
end

function Gini_fn(gg, c, dab_aux)
    S_c  = cumsum(gg .* dab_aux .* c[:])/sum(gg .* dab_aux .* c[:])
    trapez_c = 0.5*(S_c[1]*gg[1]*dab_aux[1] +
                      sum((S_c[2:end] + S_c[1:end-1]).*gg[2:end].*dab_aux[2:end]))
    C_Gini = 1 - 2*trapez_c
    return C_Gini
end

"""
```
function gstate(G1, impact; pick = Matrix{Float64}(undef, 0, 0))
```
#function  [GS,pickn,ok,uis,vs]=gstate(G1,impact,pick)
# Arguments
- `G1`:     the coefficient on  lagged y from the output of gensys
- `impact`: the coefficient on exogenous shocks from the output of gensys
- `pick`:   an optional guess at a matrix of coefficients that extracts the state vector from y

# Output
- `ok`: 0 if pick matrix is just no good.
    1 if pick matrix is usable forward, but loses initial date information
    2 if pick matrix is not usable forward, is part of a correct state vector, but is not complete
    3 if pick matrix is     usable forward, is part of a correct state vector, but is not complete
    4 if pick matrix is not usable forward, summarizes past, but is redundant
    5 if pick matrix is     usable forward, summarizes past, but is redundant
    6 if pick matrix summarizes past and is not redundant, but is not usable forward
    7 if pick matrix is perfect, both forward and as history summary

- `pickn`: a matrix of coefficients that extracts the state vector from y.  Equal
    to pick if pick is supplied and ok=1; otherwise, an ok-maximizing pick matrix that
    is somewhat similar to pick.

- `GS`: the matrix of coefficients on the lagged state variable
- `uis,vs`: If uis'*vs is square and full rank, ok=7 is possible, otherwise not.
    If ok<2, an ok>2 can be found by trying a pick in the row space of vs'.
    Any pick with pick*uis full column rank will provide a foward state (i.e. ok an odd number).
- `uiu`:   uiu'y(t)=0, and this is interpretable as the decision rule setting controls
    as functions of the states.

The solution was in the form y(t)=G1*y(t-1)+impact*z(t).
Now it's in the form y(t)=GS*pickn*y(t-1)+impact*z(t).
In general, `pickn` is mxn with m<n.
"""
function gstate(G1, impact; pick = Matrix{Float64}(undef, 0, 0))
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
    ok = oknow + 2*pinv + 4*vinp
    return GS, pickn, ok, uis, uiu, vs
end
