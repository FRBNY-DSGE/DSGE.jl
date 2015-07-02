using Base.Test

using DSGE
using DSGE: M990
include("util.jl")

# src/solve/solve.jl
function test_solve()
    model = Model()
    G1, C, impact, fmat, fwt, ywt, gev, eu, loose = solve(model)

    # Check output matrices against Matlab output (ε = 1e-4)
    G1_matlab = readcsv("models/m990/gensys/G1.csv")
    println("### G1")
    @test test_matrix_eq(G1_matlab, G1)

    C_matlab = readcsv("models/m990/gensys/C.csv")
    println("### C")
    @test test_matrix_eq(C_matlab, C)

    impact_matlab = readcsv("models/m990/gensys/impact.csv")
    println("### impact")
    @test test_matrix_eq(impact_matlab, impact)

    fmat_matlab = readcsv_complex("models/m990/gensys/fmat.csv")
    println("### fmat")
    @test test_matrix_eq(fmat_matlab, fmat)

    fwt_matlab = readcsv_complex("models/m990/gensys/fwt.csv")
    println("### fwt")
    @test test_matrix_eq(fwt_matlab, fwt)

    ywt_matlab = readcsv_complex("models/m990/gensys/ywt.csv")
    println("### ywt")
    @test test_matrix_eq(ywt_matlab, ywt)

    gev_matlab = readcsv_complex("models/m990/gensys/gev.csv")
    println("### gev")
    @test test_matrix_eq(gev_matlab, gev)

    eu_matlab = readcsv("models/m990/gensys/eu.csv")
    eu_matlab = convert(Array{Int64, 1}, vec(eu_matlab))
    println("### eu")
    @test eu_matlab == eu
end
