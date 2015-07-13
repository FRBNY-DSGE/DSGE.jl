using DSGE.Gensys



### CALLS QZDIV
function gensys_qzdiv(F::Base.LinAlg.GeneralizedSchur, C, Ψ, Π, div)
    F = schurfact(Γ0, Γ1)
    C, Ψ, Π = copy(C), copy(Ψ), copy(Π)
    
    eu = [0, 0]
    ε = 1e-6  # small number to check convergence
    nunstab = 0.0
    zxz = 0
    a, b, q, z = F[:S], F[:T], F[:Q]', F[:Z]
    n = size(a, 1)

    for i=1:n
        nunstab += (abs(b[i, i]) > div * abs(a[i,i]))
        if abs(a[i, i]) < ε && abs(b[i, i]) < ε
            zxz = 1
        end
    end

    if zxz == 1
        msg = "Coincident zeros.  Indeterminacy and/or nonexistence."
        #throw(InexactError(msg))
        throw(InexactError())
    end

    # TODO: Replicate qzdiv functionality (which we don't have permission from Chris Sims to
    #   use) using ordschur. Currently line 221 `Γ0I = inv(Γ0)` throws an exception.
    #select = abs(F[:values]) .< div
    #FS = ordschur(F, select)
    #a, b, q, z = FS[:S], FS[:T], FS[:Q]', FS[:Z]
    a, b, q, z = Gensys.qzdiv!(div, a, b, q, z)
    gev = [diag(a) diag(b)]

    q1 = q[1:n-nunstab, :]
    q2 = q[n-nunstab+1:n, :]
    z1 = z[:, 1:n-nunstab]'
    z2 = z[:, n-nunstab+1:n]'
    a2 = a[n-nunstab+1:n, n-nunstab+1:n]
    b2 = b[n-nunstab+1:n, n-nunstab+1:n]
    etawt = q2 * Π
    neta = size(Π, 2)

    # branch below is to handle case of no stable roots, which previously
    # (5/9/09) quit with an error in that case.
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
        bigev = find(diag(deta[1:md,1:md]) .> ε)
        ueta = ueta[:, bigev]
        veta = veta[:, bigev]
        deta = deta[bigev, bigev]
    end

    eu[1] = length(bigev) >= nunstab

    # eu[1] == 1 && info("gensys: Existence of a solution!")


    # ----------------------------------------------------
    # Note that existence and uniqueness are not just matters of comparing
    # numbers of roots and numbers of endogenous errors.  These counts are
    # reported below because usually they point to the source of the problem.
    # ------------------------------------------------------

    # branch below to handle case of no stable roots
    if nunstab == n
        etawt1 = zeros(0, neta)
        bigev = 0
        ueta1 = zeros(0, 0)
        veta1 = zeros(neta, 0)
        deta1 = zeros(0, 0)
    else
        etawt1 = q1 * Π
        ndeta1 = min(n-nunstab, neta)
        ueta1, deta1, veta1 = svd(etawt1)
        deta1 = diagm(deta1)  # TODO: do we need to do this
        md = min(size(deta1)...)
        bigev = find(diag(deta1[1:md, 1:md]) .> ε)
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
        nloose = sum(abs(diag(dl)) .> ε*n)
        unique = (nloose == 0)
    end

    if unique
        # info("gensys: Unique solution!")
        eu[2] = 1
    # else
    #     println("Indeterminacy. $nloose loose endog errors.")
    end

    # Cast as an int because we use it as an int!
    nunstab = int(nunstab)
    tmat = [eye(int(n-nunstab)) -(ueta*(deta\veta')*veta1*deta1*ueta1')']
    Γ0 = [tmat*a; zeros(nunstab,n-nunstab) eye(nunstab)]
    Γ1 = [tmat*b; zeros(nunstab,n)]

    # ----------------------
    # Γ0 is always non-singular because by construction there are no zeros on
    # the diagonal of a(1:n-nunstab,1:n-nunstab), which forms Γ0's ul corner.
    # -----------------------
    Γ0I = inv(Γ0)
    Γ1 = Γ0I*Γ1
    usix = n-nunstab+1:n
    C = Γ0I * [tmat*q*C; (a[usix, usix] - b[usix,usix])\q2*C]
    impact = Γ0I * [tmat*q*Ψ; zeros(nunstab, size(Ψ, 2))]
    fmat = b[usix, usix]\a[usix,usix]
    fwt = -b[usix, usix]\q2*Ψ
    ywt = Γ0I[:, usix]

    loose = Γ0I * [etawt1 * (eye(neta) - veta * veta'); zeros(nunstab, neta)]

    # -------------------- above are output for system in terms of z'y -------
    Γ1 = real(z*Γ1*z')
    C = real(z*C)
    impact = real(z * impact)
    loose = real(z * loose)

    ywt=z*ywt

    return Γ1, C, impact, fmat, fwt, ywt, gev, eu, loose
end



### CALLS ORDSCHUR
function gensys_ordschur(F::Base.LinAlg.GeneralizedSchur, C, Ψ, Π, div)
    F = schurfact(Γ0, Γ1)
    C, Ψ, Π = copy(C), copy(Ψ), copy(Π)
    
    eu = [0, 0]
    ε = 1e-6  # small number to check convergence
    nunstab = 0.0
    zxz = 0
    a, b, q, z = F[:S], F[:T], F[:Q]', F[:Z]
    n = size(a, 1)

    for i=1:n
        nunstab += (abs(b[i, i]) > div * abs(a[i,i]))
        if abs(a[i, i]) < ε && abs(b[i, i]) < ε
            zxz = 1
        end
    end

    if zxz == 1
        msg = "Coincident zeros.  Indeterminacy and/or nonexistence."
        #throw(InexactError(msg))
        throw(InexactError())
    end

    # TODO: Replicate qzdiv functionality (which we don't have permission from Chris Sims to
    #   use) using ordschur. Currently line 221 `Γ0I = inv(Γ0)` throws an exception.
    select = abs(F[:values]) .< div
    FS = ordschur(F, select)
    a, b, q, z = FS[:S], FS[:T], FS[:Q]', FS[:Z]
    #a, b, q, z = Gensys.qzdiv!(div, a, b, q, z)
    gev = [diag(a) diag(b)]

    q1 = q[1:n-nunstab, :]
    q2 = q[n-nunstab+1:n, :]
    z1 = z[:, 1:n-nunstab]'
    z2 = z[:, n-nunstab+1:n]'
    a2 = a[n-nunstab+1:n, n-nunstab+1:n]
    b2 = b[n-nunstab+1:n, n-nunstab+1:n]
    etawt = q2 * Π
    neta = size(Π, 2)

    # branch below is to handle case of no stable roots, which previously
    # (5/9/09) quit with an error in that case.
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
        bigev = find(diag(deta[1:md,1:md]) .> ε)
        ueta = ueta[:, bigev]
        veta = veta[:, bigev]
        deta = deta[bigev, bigev]
    end

    eu[1] = length(bigev) >= nunstab

    # eu[1] == 1 && info("gensys: Existence of a solution!")


    # ----------------------------------------------------
    # Note that existence and uniqueness are not just matters of comparing
    # numbers of roots and numbers of endogenous errors.  These counts are
    # reported below because usually they point to the source of the problem.
    # ------------------------------------------------------

    # branch below to handle case of no stable roots
    if nunstab == n
        etawt1 = zeros(0, neta)
        bigev = 0
        ueta1 = zeros(0, 0)
        veta1 = zeros(neta, 0)
        deta1 = zeros(0, 0)
    else
        etawt1 = q1 * Π
        ndeta1 = min(n-nunstab, neta)
        ueta1, deta1, veta1 = svd(etawt1)
        deta1 = diagm(deta1)  # TODO: do we need to do this
        md = min(size(deta1)...)
        bigev = find(diag(deta1[1:md, 1:md]) .> ε)
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
        nloose = sum(abs(diag(dl)) .> ε*n)
        unique = (nloose == 0)
    end

    if unique
        # info("gensys: Unique solution!")
        eu[2] = 1
    # else
    #     println("Indeterminacy. $nloose loose endog errors.")
    end

    # Cast as an int because we use it as an int!
    nunstab = int(nunstab)
    tmat = [eye(int(n-nunstab)) -(ueta*(deta\veta')*veta1*deta1*ueta1')']
    Γ0 = [tmat*a; zeros(nunstab,n-nunstab) eye(nunstab)]
    Γ1 = [tmat*b; zeros(nunstab,n)]

    # ----------------------
    # Γ0 is always non-singular because by construction there are no zeros on
    # the diagonal of a(1:n-nunstab,1:n-nunstab), which forms Γ0's ul corner.
    # -----------------------
    Γ0I = inv(Γ0)
    Γ1 = Γ0I*Γ1
    usix = n-nunstab+1:n
    C = Γ0I * [tmat*q*C; (a[usix, usix] - b[usix,usix])\q2*C]
    impact = Γ0I * [tmat*q*Ψ; zeros(nunstab, size(Ψ, 2))]
    fmat = b[usix, usix]\a[usix,usix]
    fwt = -b[usix, usix]\q2*Ψ
    ywt = Γ0I[:, usix]

    loose = Γ0I * [etawt1 * (eye(neta) - veta * veta'); zeros(nunstab, neta)]

    # -------------------- above are output for system in terms of z'y -------
    Γ1 = real(z*Γ1*z')
    C = real(z*C)
    impact = real(z * impact)
    loose = real(z * loose)

    ywt=z*ywt

    return Γ1, C, impact, fmat, fwt, ywt, gev, eu, loose
end
