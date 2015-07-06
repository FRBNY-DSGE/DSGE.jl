# bring things from around `Base` that we use below
using Base.LinAlg: BlasFloat, BlasInt, chkstride1, chksquare, GeneralizedSchur, Schur, blas_int
using Base.LinAlg.LAPACK: liblapack, LAPACKException
using Base: blasfunc

# The code below is pull directly from lapack.jl and schur.jl in src/linalg
# in Base. It is here to provide the ordered qz functionality present in julia
# 0.4-dev to the latest 0.3 release.

# I needed to copy/paste these here because if I just import them from Base
# they point to the function `info` instead of the local variable
macro assertargsok() #Handle only negative info codes - use only if positive info code is useful!
    :(info[1]<0 && throw(ArgumentError("invalid argument #$(-info[1]) to LAPACK call")))
end
macro lapackerror() #Handle all nonzero info codes
    :(info[1]>0 ? throw(LAPACKException(info[1])) : @assertargsok )
end

if VERSION <= v"0.4-"
    # Reorder Schur forms
    for (trsen, tgsen, elty) in
        ((:dtrsen_, :dtgsen_, :Float64),
         (:strsen_, :stgsen_, :Float32))
        @eval begin
            function trsen!(select::StridedVector{BlasInt}, T::StridedMatrix{$elty}, Q::StridedMatrix{$elty})
    # *     .. Scalar Arguments ..
    #       CHARACTER          COMPQ, JOB
    #       INTEGER            INFO, LDQ, LDT, LIWORK, LWORK, M, N
    #       DOUBLE PRECISION   S, SEP
    # *     ..
    # *     .. Array Arguments ..
    #       LOGICAL            SELECT( * )
    #       INTEGER            IWORK( * )
    #       DOUBLE PRECISION   Q( LDQ, * ), T( LDT, * ), WI( * ), WORK( * ), WR( * )
                chkstride1(T, Q)
                n = chksquare(T)
                ldt = max(1, stride(T, 2))
                ldq = max(1, stride(Q, 2))
                wr = similar(T, $elty, n)
                wi = similar(T, $elty, n)
                m = sum(select)
                work = Array($elty, 1)
                lwork = blas_int(-1)
                iwork = Array(BlasInt, 1)
                liwork = blas_int(-1)
                info = Array(BlasInt, 1)
                select = convert(Array{BlasInt}, select)

                for i = 1:2
                    ccall(($(blasfunc(trsen)), liblapack), Void,
                        (Ptr{UInt8}, Ptr{UInt8}, Ptr{BlasInt}, Ptr{BlasInt},
                        Ptr{$elty}, Ptr{BlasInt}, Ptr{$elty}, Ptr{BlasInt},
                        Ptr{$elty}, Ptr{$elty}, Ptr{BlasInt}, Ptr{Void}, Ptr{Void},
                        Ptr{$elty}, Ptr{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt},
                        Ptr{BlasInt}),
                        &'N', &'V', select, &n,
                        T, &ldt, Q, &ldq,
                        wr, wi, &m, C_NULL, C_NULL,
                        work, &lwork, iwork, &liwork,
                        info)
                    @lapackerror
                    if i == 1 # only estimated optimal lwork, liwork
                        lwork  = blas_int(real(work[1]))
                        liwork = blas_int(real(iwork[1]))
                        work   = Array($elty, lwork)
                        iwork  = Array(BlasInt, liwork)
                    end
                end
                T, Q, all(wi .== 0) ? wr : complex(wr, wi)
            end
            function tgsen!(select::StridedVector{BlasInt}, S::StridedMatrix{$elty}, T::StridedMatrix{$elty},
                                                Q::StridedMatrix{$elty}, Z::StridedMatrix{$elty})
    # *       .. Scalar Arguments ..
    # *       LOGICAL            WANTQ, WANTZ
    # *       INTEGER            IJOB, INFO, LDA, LDB, LDQ, LDZ, LIWORK, LWORK,
    # *      $                   M, N
    # *       DOUBLE PRECISION   PL, PR
    # *       ..
    # *       .. Array Arguments ..
    # *       LOGICAL            SELECT( * )
    # *       INTEGER            IWORK( * )
    # *       DOUBLE PRECISION   A( LDA, * ), ALPHAI( * ), ALPHAR( * ),
    # *      $                   B( LDB, * ), BETA( * ), DIF( * ), Q( LDQ, * ),
    # *      $                   WORK( * ), Z( LDZ, * )
    # *       ..
                chkstride1(S, T, Q, Z)
                n, nt, nq, nz = chksquare(S, T, Q, Z)
                if n != nt
                    throw(DimensionMismatch("Dimensions of S, ($n,$n), and T, ($nt,$nt), must match"))
                end
                if n != nq
                    throw(DimensionMismatch("Dimensions of S, ($n,$n), and Q, ($nq,$nq), must match"))
                end
                if n != nz
                    throw(DimensionMismatch("Dimensions of S, ($n,$n), and Z, ($nz,$nz), must match"))
                end
                lds = max(1, stride(S, 2))
                ldt = max(1, stride(T, 2))
                ldq = max(1, stride(Q, 2))
                ldz = max(1, stride(Z, 2))
                m = sum(select)
                alphai = similar(T, $elty, n)
                alphar = similar(T, $elty, n)
                beta = similar(T, $elty, n)
                lwork = blas_int(-1)
                work = Array($elty, 1)
                liwork = blas_int(-1)
                iwork = Array(BlasInt, 1)
                info = Array(BlasInt, 1)
                select = convert(Array{BlasInt}, select)

                for i = 1:2
                    ccall(($(blasfunc(tgsen)), liblapack), Void,
                           (Ptr{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt},
                            Ptr{BlasInt}, Ptr{$elty}, Ptr{BlasInt}, Ptr{$elty},
                            Ptr{BlasInt}, Ptr{$elty}, Ptr{$elty}, Ptr{$elty},
                            Ptr{$elty}, Ptr{BlasInt}, Ptr{$elty}, Ptr{BlasInt},
                            Ptr{BlasInt}, Ptr{Void}, Ptr{Void}, Ptr{Void},
                            Ptr{$elty}, Ptr{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt},
                            Ptr{BlasInt}),
                        &0, &1, &1, select,
                        &n, S, &lds, T,
                        &ldt, alphar, alphai, beta,
                        Q, &ldq, Z, &ldz,
                        &m, C_NULL, C_NULL, C_NULL,
                        work, &lwork, iwork, &liwork,
                        info)
                    @lapackerror
                    if i == 1 # only estimated optimal lwork, liwork
                        lwork  = blas_int(real(work[1]))
                        work   = Array($elty, lwork)
                        liwork = blas_int(real(iwork[1]))
                        iwork = Array(BlasInt, liwork)
                    end
                end
                S, T, complex(alphar, alphai), beta, Q, Z
            end
        end
    end

    for (trsen, tgsen, elty) in
        ((:ztrsen_, :ztgsen_, :Complex128),
         (:ctrsen_, :ctgsen_, :Complex64))
        @eval begin
            function trsen!(select::StridedVector{BlasInt}, T::StridedMatrix{$elty}, Q::StridedMatrix{$elty})
    # *     .. Scalar Arguments ..
    #       CHARACTER          COMPQ, JOB
    #       INTEGER            INFO, LDQ, LDT, LWORK, M, N
    #       DOUBLE PRECISION   S, SEP
    # *     ..
    # *     .. Array Arguments ..
    #       LOGICAL            SELECT( * )
    #       COMPLEX            Q( LDQ, * ), T( LDT, * ), W( * ), WORK( * )
                chkstride1(T, Q)
                n = chksquare(T)
                ldt = max(1, stride(T, 2))
                ldq = max(1, stride(Q, 2))
                w = similar(T, $elty, n)
                m = sum(select)
                work = Array($elty, 1)
                lwork = blas_int(-1)
                info = Array(BlasInt, 1)
                select = convert(Array{BlasInt}, select)

                for i = 1:2
                    ccall(($(blasfunc(trsen)), liblapack), Void,
                        (Ptr{UInt8}, Ptr{UInt8}, Ptr{BlasInt}, Ptr{BlasInt},
                        Ptr{$elty}, Ptr{BlasInt}, Ptr{$elty}, Ptr{BlasInt},
                        Ptr{$elty}, Ptr{BlasInt}, Ptr{Void}, Ptr{Void},
                        Ptr{$elty}, Ptr{BlasInt},
                        Ptr{BlasInt}),
                        &'N', &'V', select, &n,
                        T, &ldt, Q, &ldq,
                        w, &m, C_NULL, C_NULL,
                        work, &lwork,
                        info)
                    @lapackerror

                    if i == 1 # only estimated optimal lwork, liwork
                        lwork  = blas_int(real(work[1]))
                        work   = Array($elty, lwork)
                    end
                end
                T, Q, w
            end
            function tgsen!(select::StridedVector{BlasInt}, S::StridedMatrix{$elty}, T::StridedMatrix{$elty},
                                                Q::StridedMatrix{$elty}, Z::StridedMatrix{$elty})
    # *       .. Scalar Arguments ..
    # *       LOGICAL            WANTQ, WANTZ
    # *       INTEGER            IJOB, INFO, LDA, LDB, LDQ, LDZ, LIWORK, LWORK,
    # *      $                   M, N
    # *       DOUBLE PRECISION   PL, PR
    # *       ..
    # *       .. Array Arguments ..
    # *       LOGICAL            SELECT( * )
    # *       INTEGER            IWORK( * )
    # *       DOUBLE PRECISION   DIF( * )
    # *       COMPLEX*16         A( LDA, * ), ALPHA( * ), B( LDB, * ),
    # *      $                   BETA( * ), Q( LDQ, * ), WORK( * ), Z( LDZ, * )
    # *       ..
                chkstride1(S, T, Q, Z)
                n, nt, nq, nz = chksquare(S, T, Q, Z)
                if n != nt
                    throw(DimensionMismatch("Dimensions of S, ($n,$n), and T, ($nt,$nt), must match"))
                end
                if n != nq
                    throw(DimensionMismatch("Dimensions of S, ($n,$n), and Q, ($nq,$nq), must match"))
                end
                if n != nz
                    throw(DimensionMismatch("Dimensions of S, ($n,$n), and Z, ($nz,$nz), must match"))
                end
                lds = max(1, stride(S, 2))
                ldt = max(1, stride(T, 2))
                ldq = max(1, stride(Q, 2))
                ldz = max(1, stride(Z, 2))
                m = sum(select)
                alpha = similar(T, $elty, n)
                beta = similar(T, $elty, n)
                lwork = blas_int(-1)
                work = Array($elty, 1)
                liwork = blas_int(-1)
                iwork = Array(BlasInt, 1)
                info = Array(BlasInt, 1)
                select = convert(Array{BlasInt}, select)

                for i = 1:2
                    ccall(($(blasfunc(tgsen)), liblapack), Void,
                           (Ptr{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt},
                            Ptr{BlasInt}, Ptr{$elty}, Ptr{BlasInt}, Ptr{$elty},
                            Ptr{BlasInt}, Ptr{$elty}, Ptr{$elty},
                            Ptr{$elty}, Ptr{BlasInt}, Ptr{$elty}, Ptr{BlasInt},
                            Ptr{BlasInt}, Ptr{Void}, Ptr{Void}, Ptr{Void},
                            Ptr{$elty}, Ptr{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt},
                            Ptr{BlasInt}),
                        &0, &1, &1, select,
                        &n, S, &lds, T,
                        &ldt, alpha, beta,
                        Q, &ldq, Z, &ldz,
                        &m, C_NULL, C_NULL, C_NULL,
                        work, &lwork, iwork, &liwork,
                        info)
                    @lapackerror
                    if i == 1 # only estimated optimal lwork, liwork
                        lwork  = blas_int(real(work[1]))
                        work   = Array($elty, lwork)
                        liwork = blas_int(real(iwork[1]))
                        iwork = Array(BlasInt, liwork)
                    end
                end
                S, T, alpha, beta, Q, Z
            end
        end
    end

    ordschur!{Ty<:BlasFloat}(Q::StridedMatrix{Ty}, T::StridedMatrix{Ty}, select::Union(Vector{Bool},BitVector)) = Schur(trsen!(convert(Vector{BlasInt}, select), T , Q)...)
    ordschur{Ty<:BlasFloat}(Q::StridedMatrix{Ty}, T::StridedMatrix{Ty}, select::Union(Vector{Bool},BitVector)) = ordschur!(copy(Q), copy(T), select)
    ordschur!{Ty<:BlasFloat}(schur::Schur{Ty}, select::Union(Vector{Bool},BitVector)) = (res=ordschur!(schur.Z, schur.T, select); schur[:values][:]=res[:values]; res)
    ordschur{Ty<:BlasFloat}(schur::Schur{Ty}, select::Union(Vector{Bool},BitVector)) = ordschur(schur.Z, schur.T, select)

    ordschur!{Ty<:BlasFloat}(S::StridedMatrix{Ty}, T::StridedMatrix{Ty}, Q::StridedMatrix{Ty}, Z::StridedMatrix{Ty}, select::Union(Vector{Bool},BitVector)) = GeneralizedSchur(tgsen!(convert(Vector{BlasInt}, select), S, T, Q, Z)...)
    ordschur{Ty<:BlasFloat}(S::StridedMatrix{Ty}, T::StridedMatrix{Ty}, Q::StridedMatrix{Ty}, Z::StridedMatrix{Ty}, select::Union(Vector{Bool},BitVector)) = ordschur!(copy(S), copy(T), copy(Q), copy(Z), select)
    ordschur!{Ty<:BlasFloat}(gschur::GeneralizedSchur{Ty}, select::Union(Vector{Bool},BitVector)) = (res=ordschur!(gschur.S, gschur.T, gschur.Q, gschur.Z, select); gschur[:alpha][:]=res[:alpha]; gschur[:beta][:]=res[:beta]; res)
    ordschur{Ty<:BlasFloat}(gschur::GeneralizedSchur{Ty}, select::Union(Vector{Bool},BitVector)) = ordschur(gschur.S, gschur.T, gschur.Q, gschur.Z, select)
end
