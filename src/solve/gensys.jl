# This code is based on a routine originally copyright Chris Sims.
# See http://sims.princeton.edu/yftp/gensys/

"""
```
gensys(Γ0, Γ1, c, ψ, π)
gensys(Γ0, Γ1, c, ψ, π, div)
gensys(F::Base.LinAlg.GeneralizedSchur, c, ψ, π, div)
```

Generate state-space solution to canonical-form DSGE model.

System given as
```
Γ0*y(t) = Γ1*y(t-1) + c + Ψ*z(t) + Π*η(t),
```
with z an exogenous variable process and η being endogenously
determined one-step-ahead expectational errors.

Returned system is
```
y(t) = G1*y(t-1) + C + impact*z(t) + ywt*inv(I-fmat*inv(L))*fwt*z(t+1)
```
Returned values are
```
G1, C, impact, fmat, fwt, ywt, gev, eu, loose
```

If `z(t)` is i.i.d., the last term drops out.

If `div` is omitted from argument list, a `div`>1 is calculated.

### Return codes

* `eu[1]==1` for existence
* `eu[2]==1` for uniqueness
* `eu[1]==-1` for existence only with not-s.c. z;
* `eu==[-2,-2]` for coincident zeros.

### Notes

We constrain Julia to use the complex version of the `schurfact` routine regardless of the
types of `Γ0` and `Γ1`, to match the behavior of Matlab.  Matlab always uses the complex version
of the Schur decomposition, even if the inputs are real numbers.
"""
function gensys(Γ0, Γ1, c, ψ, π)
    F = schurfact(complex(Γ0), complex(Γ1))
    div = new_div(F)
    gensys(F, c, ψ, π, div)
end
function gensys(Γ0, Γ1, c, ψ, π, div)
    F = schurfact(complex(Γ0), complex(Γ1))
    gensys(F, c, ψ, π, div)
end

# Method that does the real work. Work directly on the decomposition F
function gensys(F::Base.LinAlg.GeneralizedSchur, c, ψ, π, div)
    eu = [0, 0]
    ϵ = 1e-6  # small number to check convergence
    nunstab = 0.0
    zxz = 0
    a, b, q, z = F[:S], F[:T], F[:Q]', F[:Z]
    n = size(a, 1)

    for i=1:n
        nunstab += (abs(b[i, i]) > div * abs(a[i,i]))
        if abs(a[i, i]) < ϵ && abs(b[i, i]) < ϵ
            zxz = 1
        end
    end

    if zxz == 1
        throw(GensysError("Coincident zeros.  Indeterminacy and/or nonexistence."))
    end

    select = (abs(F[:alpha]) .> div*abs(F[:beta]))
    FS = ordschur(F, select)
    a, b, q, z = FS[:S], FS[:T], FS[:Q]', FS[:Z]
    gev = [diag(a) diag(b)]

    nunstab_int = round(Int,nunstab)

    q1 = q[1:n-nunstab_int, :]
    q2 = q[n-nunstab_int+1:n, :]
    z1 = z[:, 1:n-nunstab_int]'
    z2 = z[:, n-nunstab_int+1:n]'
    a2 = a[n-nunstab_int+1:n, n-nunstab_int+1:n]
    b2 = b[n-nunstab_int+1:n, n-nunstab_int+1:n]

    etawt = q2 * π
    neta = size(π, 2)

    # branch below is to handle case of no stable roots, rather than quitting with an error
    # in that case.
    if isapprox(nunstab, 0.0)
        etawt == zeros(0, neta)
        ueta = zeros(0, 0)
        deta = zeros(0, 0)
        veta = zeros(neta, 0)
        bigev = 0
    else
        ueta, deta, veta = svd(etawt)
        deta = diagm(deta)  # TODO: do we need to do this
        md = min(size(deta)...)
        bigev = find(diag(deta[1:md,1:md]) .> ϵ)
        ueta = ueta[:, bigev]
        veta = veta[:, bigev]
        deta = deta[bigev, bigev]
    end

    eu[1] = length(bigev) >= nunstab

    # eu[1] == 1 && info("gensys: Existence of a solution!")

    # Note that existence and uniqueness are not just matters of comparing
    # numbers of roots and numbers of endogenous errors.  These counts are
    # reported below because usually they point to the source of the problem.

    # branch below to handle case of no stable roots
    if nunstab == n
        etawt1 = zeros(0, neta)
        bigev = 0
        ueta1 = zeros(0, 0)
        veta1 = zeros(neta, 0)
        deta1 = zeros(0, 0)
    else
        etawt1 = q1 * π
        ndeta1 = min(n-nunstab, neta)
        ueta1, deta1, veta1 = svd(etawt1)
        deta1 = diagm(deta1)  # TODO: do we need to do this
        md = min(size(deta1)...)
        bigev = find(diag(deta1[1:md, 1:md]) .> ϵ)
        ueta1 = ueta1[:, bigev]
        veta1 = veta1[:, bigev]
        deta1 = deta1[bigev, bigev]
    end

    if isempty(veta1)
        unique = 1
    else
        loose = veta1-veta*veta'*veta1
        ul, dl, vl = svd(loose)
        dl = diagm(dl)  # TODO: do we need to do this
        nloose = sum(abs(diag(dl)) .> ϵ*n)
        unique = (nloose == 0)
    end

    if unique
        # info("gensys: Unique solution!")
        eu[2] = 1
    # else
    #     println("Indeterminacy. $nloose loose endog errors.")
    end

    # Cast as an int because we use it as an int!
    nunstab = round(Int, nunstab)
    tmat = [eye(convert(Int64,(n-nunstab))) -(ueta*(deta\veta')*veta1*deta1*ueta1')']

    G0 = [tmat*a; zeros(nunstab,n-nunstab) eye(nunstab)]
    G1 = [tmat*b; zeros(nunstab,n)]

    # G0 is always non-singular because by construction there are no zeros on
    # the diagonal of a(1:n-nunstab,1:n-nunstab), which forms G0's ul corner.
    G0I = inv(G0)
    G1 = G0I*G1
    usix = n-nunstab+1:n
    C = G0I * [tmat*q*c; (a[usix, usix] - b[usix,usix])\q2*c]
    impact = G0I * [tmat*q*ψ; zeros(nunstab, size(ψ, 2))]
    fmat = b[usix, usix]\a[usix,usix]
    fwt = -b[usix, usix]\q2*ψ
    ywt = G0I[:, usix]

    loose = G0I * [etawt1 * (eye(neta) - veta * veta'); zeros(nunstab, neta)]

    # above are output for system in terms of z'y
    G1 = real(z*G1*z')
    C = real(z*C)
    impact = real(z * impact)
    loose = real(z * loose)

    ywt=z*ywt
   
    #Commented out since Schorfheide's code allows for indeterminacy/non-uniqueness
    #if eu[1] != 1 || eu[2] != 1
    #    throw(GensysError("Gensys does not give existence and uniqueness."))
    #end

    return G1, C, impact, fmat, fwt, ywt, gev, eu, loose
end

"""
```
GensysError <: Exception
```

A `GensysError` is thrown when Gensys does not give a unique solution, or no solution
exists. If a `GensysError`is thrown during Metropolis-Hastings, it is caught by `posterior`.
`posterior` then returns a value of `-Inf`, which Metropolis-Hastings always rejects.

### Fields

* `msg::ASCIIString`: Info message. Default = "Error in gensys." 
"""
type GensysError <: Exception
    msg::ASCIIString
end
GensysError() = GensysError("Error in gensys.")
Base.showerror(io::IO, ex::GensysError) = print(io, ex.msg)

function new_div(F::Base.LinAlg.GeneralizedSchur)
    ϵ = 1e-6  # small number to check convergence
    n = size(F[:T], 1)

    a, b, q, z = F[:S], F[:T], F[:Q]', F[:Z]

    div = 1.01

    for i=1:n
        if abs(a[i, i]) > 0
            divhat = abs(b[i, i]) / abs(a[i, i])
            if 1 + ϵ < divhat && divhat <= div
                div = .5 * (1 + divhat)
            end
        end
    end

    return div
end
