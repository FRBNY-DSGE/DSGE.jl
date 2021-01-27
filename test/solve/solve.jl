using HDF5, Test, DSGE, ModelConstructors

path = dirname(@__FILE__)

file = "$path/../reference/solve.h5"
TTT_expected = h5read(file, "TTT")
CCC_expected = h5read(file, "CCC")
RRR_expected = h5read(file, "RRR")

m = AnSchorfheide()
TTT, RRR, CCC = solve(m)

@testset "Check state-space system matches reference" begin
    @test @test_matrix_approx_eq TTT_expected TTT
    @test @test_matrix_approx_eq RRR_expected RRR
    @test @test_matrix_approx_eq CCC_expected CCC
end

m = Model1002("ss10")
m <= Setting(:regime_switching, true)
setup_permanent_altpol!(m, taylor93())
TTT_alt93, RRR_alt93, CCC_alt93 = solve(m; regime_switching = true, regimes = collect(1:2), gensys_regimes = [1:2])
file = "$path/../reference/solve_alt.h5"
TTT_alt_expected93 = h5read(file, "TTT_alt93")
CCC_alt_expected93 = h5read(file, "CCC_alt93")
RRR_alt_expected93 = h5read(file, "RRR_alt93")

@testset "Check state-space system with altpolicy matches reference for Taylor93" begin
    @test @test_matrix_approx_eq TTT_alt_expected93 TTT_alt93[2]
    @test @test_matrix_approx_eq RRR_alt_expected93 RRR_alt93[2]
    @test @test_matrix_approx_eq CCC_alt_expected93 CCC_alt93[2]
    @test maximum(abs.(TTT_alt_expected93 - TTT_alt93[1])) > 1e-6
    @test maximum(abs.(RRR_alt_expected93 - RRR_alt93[1])) > 1e-6
end

m = Model1002("ss10")
m <= Setting(:regime_switching, true)
setup_permanent_altpol!(m, taylor99())
TTT_alt99, RRR_alt99, CCC_alt99 = solve(m; regime_switching = true, regimes = collect(1:2), gensys_regimes = [1:2])
file = "$path/../reference/solve_alt.h5"
TTT_alt_expected99 = h5read(file, "TTT_alt99")
CCC_alt_expected99 = h5read(file, "CCC_alt99")
RRR_alt_expected99 = h5read(file, "RRR_alt99")

@testset "Check state-space system with altpolicy matches reference for Taylor99" begin
    @test @test_matrix_approx_eq TTT_alt_expected99 TTT_alt99[2]
    @test @test_matrix_approx_eq RRR_alt_expected99 RRR_alt99[2]
    @test @test_matrix_approx_eq CCC_alt_expected99 CCC_alt99[2]
    @test maximum(abs.(TTT_alt_expected99 - TTT_alt99[1])) > 1e-6
    @test maximum(abs.(RRR_alt_expected99 - RRR_alt99[1])) > 1e-6
end

nothing
