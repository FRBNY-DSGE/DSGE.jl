include("../src/ordered_qz.jl")

debug = false

using Base.Test

# these tests were pulled directly from main Julia
import Base.LinAlg: BlasComplex, BlasFloat, BlasReal, QRPivoted

n = 10

# Split n into 2 parts for tests needing two matrices
n1 = div(n, 2)
n2 = 2*n1

srand(1234321)
a = rand(n,n)
areal = randn(n,n)/2
aimg  = randn(n,n)/2
a2real = randn(n,n)/2
a2img  = randn(n,n)/2
breal = randn(n,2)/2
bimg  = randn(n,2)/2


for eltya in (Float32, Float64, Complex64, Complex128, BigFloat, Int)
    a = eltya == Int ? rand(1:7, n, n) : convert(Matrix{eltya}, eltya <: Complex ? complex(areal, aimg) : areal)
    a2 = eltya == Int ? rand(1:7, n, n) : convert(Matrix{eltya}, eltya <: Complex ? complex(a2real, a2img) : a2real)
    asym = a'+a                  # symmetric indefinite
    apd  = a'*a                 # symmetric positive-definite
    ε = εa = eps(abs(float(one(eltya))))

    for eltyb in (Float32, Float64, Complex64, Complex128, Int)
        b = eltyb == Int ? rand(1:5, n, 2) : convert(Matrix{eltyb}, eltyb <: Complex ? complex(breal, bimg) : breal)
        εb = eps(abs(float(one(eltyb))))
        ε = max(εa,εb)

debug && println("Reorder Schur")
    if eltya != BigFloat && eltyb != BigFloat # Revisit when implemented in julia
        # use asym for real schur to enforce tridiag structure
        # avoiding partly selection of conj. eigenvalues
        ordschura = eltya <: Complex ? a : asym
        S = schurfact(ordschura)
        select = bitrand(n)
        O = ordschur(S, select)
        bool(sum(select)) && @test_approx_eq S[:values][find(select)] O[:values][1:sum(select)]
        @test_approx_eq O[:vectors]*O[:Schur]*O[:vectors]' ordschura
    end

debug && println("Reorder Generalized Schur")
    if eltya != BigFloat && eltyb != BigFloat # Revisit when implemented in Julia
        a1_sf = a[1:n1, 1:n1]
        a2_sf = a[n1+1:n2, n1+1:n2]
        NS = schurfact(a1_sf, a2_sf)
        # Currently just testing with selecting gen eig values < 1
        select = abs2(NS[:values]) .< 1
        m = sum(select)
        S = ordschur(NS, select)
        # Make sure that the new factorization stil factors matrix
        @test_approx_eq S[:Q]*S[:S]*S[:Z]' a1_sf
        @test_approx_eq S[:Q]*S[:T]*S[:Z]' a2_sf
        # Make sure that we have sorted it correctly
        @test_approx_eq NS[:values][find(select)] S[:values][1:m]
    end

end  # eltyb
end  # eltya
