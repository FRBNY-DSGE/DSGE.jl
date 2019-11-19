"""
```
gensysct(Γ0, Γ1, c, Ψ, Π; ϵ, div)
gensysct(F::LinearAlgebra.GeneralizedSchur, c, Ψ, Π; ϵ, div)
```

Generate state-space solution to canonical-form DSGE model.

System given as
```
Γ0*Dy(t) = Γ1*y(t) + c + Ψ*z(t) + Π*η(t),
```
with z an exogenous variable process and η being endogenously white noise expectational errors.

Returned system is
```
Dy(t) = G1*y(t) + C + impact * z(t)
```


Returned values are
```
G1, C, impact, qt', a, b, z, eu
```

This code is based on a [routine](http://sims.princeton.edu/yftp/gensys/) originally copyright Chris Sims.

Here for completeness, but should not be used for most HACT/HANK models.

Also returned is the qz decomposition, qt'az' = Γ0, qt'bz' = Γ1, with a and b
upper triangular and the system ordered so that all zeros on the diagonal of b are in
the lower right corner, all cases where the real part of bii/aii is greater than or
equal to div appear in the next block above the zeros, and the remaining bii/aii's
all have bii/aii<div .  These elements can be used to construct the full backward and
forward solution.  See the paper \"Solving Linear Rational Expectations Models\",
http://eco-072399b.princeton.edu/yftp/gensys .  Note that if one simply wants the backwar
and forward projection of y on eΨlon, ignoring existence and uniqueness questions, the
projection can be computed by Fourier methods.

If `div` is omitted from argument list, a `div`>1 is calculated.

### Return codes

* `eu[1]==1` for existence
* `eu[2]==1` for uniqueness
* `eu[1]==-1` for existence for white noise η
* `eu==[-2,-2]` for coincident zeros.
"""
function gensysct(Γ0::Matrix{Float64}, Γ1::Matrix{Float64}, c::Array{Float64},
                  Ψ::Matrix{Float64}, Π::Matrix{Float64};
                  check_existence::Bool = true, check_uniqueness::Bool = true,
                  ϵ::Float64 = sqrt(eps()) * 10,
                  div::Float64 = -1.0, complex_decomposition::Bool = false)

    F = complex_decomposition ? schurfact(complex(Γ0), complex(Γ1)) : schurfact(real(Γ0), real(Γ1))
    gensysct(F, c, Ψ, Π; check_existence = check_existence, check_uniqueness = check_uniqueness,
             ϵ = ϵ, div = div)
end

# In place version. Changes Γ0, Γ1
function gensysct!(Γ0::Matrix{Float64}, Γ1::Matrix{Float64}, c::Array{Float64},
                   Ψ::Matrix{Float64}, Π::Matrix{Float64};
                   check_existence::Bool = true, check_uniqueness::Bool = true,
                   ϵ::Float64 = sqrt(eps()) * 10,
                   div::Float64 = -1.0, complex_decomposition::Bool = false)

    F = complex_decomposition ? schur!(complex(Γ0), complex(Γ1)) : schur!(real(Γ0), real(Γ1))
    gensysct(F, c, Ψ, Π; check_existence = check_existence, check_uniqueness = check_uniqueness,
             ϵ = ϵ, div = div)
end

# Uses basic schur decomposition
function gensysct(Γ1::Matrix{Float64}, c::Array{Float64}, Ψ::Matrix{Float64}, Π::Matrix{Float64};
                  check_existence::Bool = true, check_uniqueness::Bool = true,
                  ϵ::Float64 = sqrt(eps()) * 10,
                  div::Float64 = -1.0, complex_decomposition::Bool = false)
    # Assumes that Γ0 is the identity
    F = complex_decomposition ? schurfact(complex(Γ1)) : schurfact(real(Γ1))
    gensysct(F, c, Ψ, Π; check_existence = check_existence, check_uniqueness = check_uniqueness,
             ϵ = ϵ, div = div)
end

# In place version of above gensysct
function gensysct!(Γ1::Matrix{Float64}, c::Array{Float64}, Ψ::Matrix{Float64}, Π::Matrix{Float64};
                  check_existence::Bool = true, check_uniqueness::Bool = true, ϵ::Float64 = sqrt(eps()) * 10,
                  div::Float64 = -1.0, complex_decomposition::Bool = false)
    # Assumes that Γ0 is the identity
    F = complex_decomposition ? schur!(complex(Γ1)) : schur!(real(Γ1))
    gensysct(F, c, Ψ, Π; check_existence = check_existence, check_uniqueness = check_uniqueness,
             ϵ = ϵ, div = div)
end

function gensysct(F::LinearAlgebra.GeneralizedSchur, c::Array{Float64}, Ψ::Matrix{Float64}, Π::Matrix{Float64};
                  check_existence::Bool = true, check_uniqueness::Bool = true, ϵ::Float64 = sqrt(eps()) * 10,
                  div::Float64 = -1.0)
    div < 0.0 && (div = new_divct(F; ϵ = ϵ))
    eu = [0, 0]
    a, b = F.S, F.T
    n = size(a, 1)

    for i in 1:n
        if (abs(a[i, i]) < ϵ) && (abs(b[i, i]) < ϵ)
            info("Coincident zeros.  Indeterminacy and/or nonexistence.")
            eu = [-2, -2]
            G1 = Array{Float64, 2}() ;  C = Array{Float64, 1}() ; impact = Array{Float64, 2}()
            a, b, qt, z = FS.S, FS.T, FS.Q, FS.Z
            return G1, C, impact, qt', a, b, z, eu
        end
    end
    movelast = fill(false, n)
    for i in 1:n
        movelast[i] = (real(b[i, i] / a[i, i]) > div) || (abs(a[i, i]) < ϵ)
    end
    nunstab = sum(movelast)
    FS = ordschur!(F, .!movelast)
    a, b, qt, z = FS.S, FS.T, FS.Q, FS.Z

    qt1 = qt[:, 1:(n - nunstab)]
    qt2 = qt[:, (n - nunstab + 1):n]
    a2 = a[(n - nunstab + 1):n, (n - nunstab + 1):n]
    b2 = b[(n - nunstab + 1):n, (n - nunstab + 1):n]
    etawt = qt2' * Π # Ac_mul_B(qt2, Π)
    ~, ueta, deta, veta = decomposition_svdct!(etawt, ϵ = ϵ)
    zwt = qt2' * Ψ # Ac_mul_B(qt2, Ψ)
    bigev, uz, dz, vz = decomposition_svdct!(zwt, ϵ = ϵ)

    # check existence
    if isempty(bigev) || !check_existence
        exist = true
        existx = true
    else
        exist = norm(uz- (ueta * adjoint(ueta)) * uz, 2) < ϵ * n

        zwtx0 = b2 \ zwt
        zwtx = zwtx0
        M = b2 \ a2
        M = scale!(M, 1 / norm(M))
        for i in 2:nunstab
            zwtx = hcat(M * zwtx, zwtx0)
        end
        zwtx = b2 * zwtx
        bigev, ux, dx, vx = decomposition_svdct!(zwtx, ϵ = ϵ)
        existx = norm(ux - ueta * ueta' * ux, 2) < ϵ * n
    end

    etawt1 = qt1' * Π # Ac_mul_B(qt1, Π)
    bigev, ueta1, deta1, veta1 = decomposition_svdct!(etawt1, ϵ = ϵ)

    if nunstab == 0
       eu[1] = 1
    else
        if exist
            eu[1] = -1
        end
    end

    # check uniqueness
    if isempty(veta1) || !check_uniqueness
        eu[2] = 1
    else
        eu[2] = Int64(norm(veta1- (veta * adjoint(veta)) * veta1, 2) < ϵ * n)
    end

    tmat = hcat(eye(n - nunstab), -ueta1 * deta1 * conj(veta1)' * veta * (deta \ ueta'))
    G0 =  vcat(tmat * a, hcat(zeros(nunstab, n - nunstab), eye(nunstab)))
    G1 =  vcat(tmat * b, zeros(nunstab, n))

    G0I = inv(G0)
    G1 = G0I * G1
    usix = (n - nunstab + 1):n
    Busix = b[usix,usix]
    Ausix = a[usix,usix]
    C = G0I * vcat(tmat * qt' * c, (Ausix - Busix) \ (qt2' * c))
    impact = G0I * vcat(tmat * qt' * Ψ, zeros(nunstab, size(Ψ, 2)))
    # C = G0I * vcat(tmat * Ac_mul_B(qt, c), (Ausix - Busix) \ Ac_mul_B(qt2, c))
    # impact = G0I * vcat(tmat * Ac_mul_B(qt, Ψ), zeros(nunstab, size(Ψ, 2)))

    #G1 = z * A_mul_Bc(G1, z)
    G1 = z * G1 * z'
    G1 = real(G1)
    C = real(z * C) .* ones(size(z, 1), 1)
    impact = real(z * impact)

    return G1::Matrix{Float64}, C::Matrix{Float64}, impact::Matrix{Float64}, qt'::Matrix, a::Matrix, b::Matrix, z::Matrix, eu::Vector{Int64}
end

# Port of SeHyoun's schur_solver
## Inputs
# div: number of unstable roots known a priori
## Outputs
# eu[1] = 1 solution exists
#         0 no stable solution
#         -1 solution exists for white noise η
# eu[2] = 1 solution is unique
#         0 solution is not unique
function gensysct(F::LinearAlgebra.Schur, c::Array{Float64}, Ψ::Matrix{Float64}, Π::Matrix{Float64};
                  check_existence::Bool = true, check_uniqueness::Bool = true, ϵ::Float64 = sqrt(eps()) * 10,
                  div::Float64 = -1.0)
    eu = [0, 0]
    U, T = F.Z, F.T        # U is unitary matrix, T is Schur matrix
    n = size(U, 1)
    g_eigs = real(F.values) # Get real eigenvalues
    stable_eigs = g_eigs .<= 0 # stable eigenvalues
    nunstab = n - sum(stable_eigs)  # Unstable eigenvalues
    if div > -1
        sort_order = sortperm(g_eigs, rev = true) # g_eigs[sort_order] -> descending sorted g_eigs
        locs = ones(n, 1)
        locs[sort_order[1:div]] = 0 # Set to zero the eigenvalues which are unstable
        FS = ordschur!(F, locs)     # Change F in place, FS creates reference to same object, so no copying
        if nunstab > div
            @warn "<gensysct>: There are more than div number of positive eigenvalues with smallest values:\n"
            println(g_eigs[sort_order][div + 1:nunstab])
        elseif nunstab < div
            @warn "<gensysct>: There are less than div number of positive eigenvalues:\n"
            println(g_eigs[sort_order][nunstab + 1:div])
        end
        nunstab = div
    else
        FS = ordschur!(F, stable_eigs)
        U, T = FS.Z, FS.T
    end

    U1 = U[:, 1:(n - nunstab)]'
    U2 = U[:, (n - nunstab + 1):n]'
    etawt = U2 * Π # Ac_mul_B(U2, Π)
    ~, ueta, deta, veta = decomposition_svdct!(etawt, ϵ = ϵ) # Wrapper function here, see end of script, this checks for bigev and returns ueta, etc., after selecting only the relevant columns

    # check existence
    if check_existence
        zwt = U2 * Ψ   # Ac_mul_B(U2, Ψ)
        bigev, uz, dz, vz = decomposition_svdct!(zwt, ϵ = ϵ)
        if isempty(bigev)
            eu[1] = 1
        else
            eu[1] = norm(uz- (ueta * adjoint(ueta)) * uz, 2) < ϵ * n
        end

        if (eu[1] == 0 && (div == -1))
            @warn "<gensysct>: Solution does not exist"
        end
        impact = real(-Π * veta * (deta \ ueta') * uz * dz * vz' + Ψ)
    else
        eu[1] = 1
        impact = real(-Π * veta * (deta \ ueta') * U2 + Ψ)
    end

    # Check uniqueness
    if check_uniqueness
        etawt1 = U1 * Π # Ac_mul_B(U1, Π)
        ~, ~, deta1, veta1 = decomposition_svdct!(etawt1)
        if isempty(veta1)
            eu[2] = 1
        else
            eu[2] = Int64(norm(veta1 - (veta * adjoint(veta)) * veta1, 2) < ϵ * n)
        end
    end
    I, J, V = SparseArrays.spdiagm_internal(0 => [ones(n - nunstab); zeros(nunstab)])
    diag_m = sparse(I, J, V, n, n)
    G1 = real(U * T * diag_m * U')#spdiagm([ones(n - nunstab); zeros(nunstab)], 0, n, n) * U')
    F = U1[:, 1:nunstab]' * inv(U1[:, nunstab + 1:end]')
    impact = [F * Ψ[nunstab + 1:end, :]; Ψ[nunstab + 1:end, :]]
    C = real(U * c) .* ones(size(U, 1), 1)

    return G1::Matrix{Float64}, C::Matrix{Float64}, impact::Matrix{Float64}, U::Matrix, T::Matrix::Matrix, eu::Vector{Int64}
end

function new_divct(F::LinearAlgebra.GeneralizedSchur; ϵ::Float64 = sqrt(eps()) * 10)
    a, b = F.S, F.T
    n = size(a, 1)
    div = 0.001
    for i in 1:n
        if abs(a[i, i]) > ϵ
            divhat = real(b[i, i] / a[i, i])
            if (ϵ < divhat) && (divhat < div)
                div = 0.5 * divhat
            end
        end
    end
    return div
end

function decomposition_svdct!(A; ϵ::Float64 = sqrt(eps()) * 10)
    Asvd = svd!(A)
    bigev = findall(Asvd.S .> ϵ)
    Au = Asvd.U[:, bigev]
    Ad = diagm(0 => Asvd.S[bigev])
    Av = Asvd.V[:, bigev]
    return bigev, Au, Ad, Av
end

# Additional cases when inputs are sparse matrices
function gensysct(Γ0::SparseMatrixCSC{Float64,Int64}, Γ1::SparseMatrixCSC{Float64,Int64},
                  c::SparseMatrixCSC{Float64,Int64}, Ψ::SparseMatrixCSC{Float64,Int64},
                  Π::SparseMatrixCSC{Float64,Int64};
                  check_existence::Bool = true, check_uniqueness::Bool = true, ϵ::Float64 = sqrt(eps()) * 10,
                  div::Float64 = -1.0, complex_decomposition::Bool = false)
    F = complex_decomposition ? schurfact(complex(Matrix{Float64}(Γ0)), complex(Matrix{Float64}(Γ1))) : schurfact(real(Matrix{Float64}(Γ0)), real(Matrix{Float64}(Γ1)))
    gensysct(F, Matrix{Float64}(c), Matrix{Float64}(Ψ), Matrix{Float64}(Π); check_existence = check_existence,
             check_uniqueness = check_uniqueness, ϵ = ϵ, div = div)
end

function gensysct!(Γ0::SparseMatrixCSC{Float64,Int64}, Γ1::SparseMatrixCSC{Float64,Int64},
                   c::SparseMatrixCSC{Float64,Int64}, Ψ::SparseMatrixCSC{Float64,Int64},
                   Π::SparseMatrixCSC{Float64,Int64};
                  check_existence::Bool = true, check_uniqueness::Bool = true, ϵ::Float64 = sqrt(eps()) * 10,
                  div::Float64 = -1.0, complex_decomposition::Bool = false)
    F = complex_decomposition ? schur!(complex(Matrix{Float64}(Γ0)), complex(Matrix{Float64}(Γ1))) : schur!(real(Matrix{Float64}(Γ0)), real(Matrix{Float64}(Γ1)))
    gensysct(F, Matrix{Float64}(c), Matrix{Float64}(Ψ), Matrix{Float64}(Π); check_existence = check_existence,
             check_uniqueness = check_uniqueness, ϵ = ϵ, div = div)
end

function gensysct(Γ1::SparseMatrixCSC{Float64,Int64}, c::SparseMatrixCSC{Float64,Int64},
                  Ψ::SparseMatrixCSC{Float64,Int64}, Π::SparseMatrixCSC{Float64,Int64};
                  check_existence::Bool = true, check_uniqueness::Bool = true, ϵ::Float64 = sqrt(eps()) * 10,
                  div::Float64 = -1.0, complex_decomposition::Bool = false)
    # Assumes that Γ0 is the identity
    F = complex_decomposition ? schurfact(complex(Matrix{Float64}(Γ1))) : schurfact(real(Matrix{Float64}(Γ1)))
    gensysct(F, Matrix{Float64}(c), Matrix{Float64}(Ψ), Matrix{Float64}(Π); check_existence = check_existence,
             check_uniqueness = check_uniqueness, ϵ = ϵ, div = div)
end

function gensysct!(Γ1::SparseMatrixCSC{Float64,Int64}, c::SparseMatrixCSC{Float64,Int64},
                   Ψ::SparseMatrixCSC{Float64,Int64}, Π::SparseMatrixCSC{Float64,Int64};
                  check_existence::Bool = true, check_uniqueness::Bool = true, ϵ::Float64 = sqrt(eps()) * 10,
                  div::Float64 = -1.0, complex_decomposition::Bool = false)
    # Assumes that Γ0 is the identity
    F = complex_decomposition ? schur!(complex(Matrix{Float64}(Γ1))) : schur!(real(Matrix{Float64}(Γ1)))
    gensysct(F, Matrix{Float64}(c), Matrix{Float64}(Ψ), Matrix{Float64}(Π); check_existence = check_existence,
             check_uniqueness = check_uniqueness, ϵ = ϵ, div = div)
end
