@testset "Check miscellaneous utility functions" begin
    @test DSGE.abbrev_symbol(:aaaaa, 3) == "aaa"
    v = [1,2,3]
    DSGE.sorted_list_insert!(v, 4)
    @test v == [1,2,3,4]
    @test DSGE.quarter_range(DSGE.quartertodate("2000-Q1"), DSGE.quartertodate("2000-Q2")) == [Date("2000-03-31",
                                                                                               "yyyy-mm-dd"),
                                                                                          Date("2000-06-30",
                                                                                               "yyyy-mm-dd")]
end

@testset "Check detexify" begin
    @test DSGE.detexify(:α) == :alpha
    @test DSGE.detexify(:β) == :beta
    @test DSGE.detexify(:γ) == :gamma
    @test DSGE.detexify(:δ) == :delta
    @test DSGE.detexify(:ϵ) == :epsilon
    @test DSGE.detexify(:ε) == :epsilon
    @test DSGE.detexify(:ζ) == :zeta
    @test DSGE.detexify(:η) == :eta
    @test DSGE.detexify(:θ) == :theta
    @test DSGE.detexify(:ι) == :iota
    @test DSGE.detexify(:κ) == :kappa
    @test DSGE.detexify(:λ) == :lambda
    @test DSGE.detexify(:μ) == :mu
    @test DSGE.detexify(:ν) == :nu
    @test DSGE.detexify(:ξ) == :xi
    @test DSGE.detexify(:π) == :pi
    @test DSGE.detexify(:ρ) == :rho
    @test DSGE.detexify(:ϱ) == :rho
    @test DSGE.detexify(:σ) == :sigma
    @test DSGE.detexify(:ς) == :sigma
    @test DSGE.detexify(:τ) == :tau
    @test DSGE.detexify(:υ) == :upsilon
    @test DSGE.detexify(:ϕ) == :phi
    @test DSGE.detexify(:φ) == :varphi
    @test DSGE.detexify(:χ) == :chi
    @test DSGE.detexify(:ψ) == :psi
    @test DSGE.detexify(:ω) == :omega

    @test DSGE.detexify(:Α) == :Alpha
    @test DSGE.detexify(:Β) == :Beta
    @test DSGE.detexify(:Γ) == :Gamma
    @test DSGE.detexify(:Δ) == :Delta
    @test DSGE.detexify(:Ε) == :Epsilon
    @test DSGE.detexify(:Ζ) == :Zeta
    @test DSGE.detexify(:Η) == :Eta
    @test DSGE.detexify(:Θ) == :Theta
    @test DSGE.detexify(:Ι) == :Iota
    @test DSGE.detexify(:Κ) == :Kappa
    @test DSGE.detexify(:Λ) == :Lambda
    @test DSGE.detexify(:Μ) == :Mu
    @test DSGE.detexify(:Ν) == :Nu
    @test DSGE.detexify(:Ξ) == :Xi
    @test DSGE.detexify(:Π) == :Pi
    @test DSGE.detexify(:Ρ) == :Rho
    @test DSGE.detexify(:Σ) == :Sigma
    @test DSGE.detexify(:Τ) == :Tau
    @test DSGE.detexify(:Υ) == :Upsilon
    @test DSGE.detexify(:Φ) == :Phi
    @test DSGE.detexify(:Χ) == :Chi
    @test DSGE.detexify(:Ψ) == :Psi
    @test DSGE.detexify(:Ω) == :Omega

    @test DSGE.detexify("α") == "alpha"
    @test DSGE.detexify("β") == "beta"
    @test DSGE.detexify("γ") == "gamma"
    @test DSGE.detexify("δ") == "delta"
    @test DSGE.detexify("ϵ") == "epsilon"
    @test DSGE.detexify("ε") == "epsilon"
    @test DSGE.detexify("ζ") == "zeta"
    @test DSGE.detexify("η") == "eta"
    @test DSGE.detexify("θ") == "theta"
    @test DSGE.detexify("ι") == "iota"
    @test DSGE.detexify("κ") == "kappa"
    @test DSGE.detexify("λ") == "lambda"
    @test DSGE.detexify("μ") == "mu"
    @test DSGE.detexify("ν") == "nu"
    @test DSGE.detexify("ξ") == "xi"
    @test DSGE.detexify("π") == "pi"
    @test DSGE.detexify("ρ") == "rho"
    @test DSGE.detexify("ϱ") == "rho"
    @test DSGE.detexify("σ") == "sigma"
    @test DSGE.detexify("ς") == "sigma"
    @test DSGE.detexify("τ") == "tau"
    @test DSGE.detexify("υ") == "upsilon"
    @test DSGE.detexify("ϕ") == "phi"
    @test DSGE.detexify("φ") == "varphi"
    @test DSGE.detexify("χ") == "chi"
    @test DSGE.detexify("ψ") == "psi"
    @test DSGE.detexify("ω") == "omega"

    @test DSGE.detexify("Α") == "Alpha"
    @test DSGE.detexify("Β") == "Beta"
    @test DSGE.detexify("Γ") == "Gamma"
    @test DSGE.detexify("Δ") == "Delta"
    @test DSGE.detexify("Ε") == "Epsilon"
    @test DSGE.detexify("Ζ") == "Zeta"
    @test DSGE.detexify("Η") == "Eta"
    @test DSGE.detexify("Θ") == "Theta"
    @test DSGE.detexify("Ι") == "Iota"
    @test DSGE.detexify("Κ") == "Kappa"
    @test DSGE.detexify("Λ") == "Lambda"
    @test DSGE.detexify("Μ") == "Mu"
    @test DSGE.detexify("Ν") == "Nu"
    @test DSGE.detexify("Ξ") == "Xi"
    @test DSGE.detexify("Π") == "Pi"
    @test DSGE.detexify("Ρ") == "Rho"
    @test DSGE.detexify("Σ") == "Sigma"
    @test DSGE.detexify("Τ") == "Tau"
    @test DSGE.detexify("Υ") == "Upsilon"
    @test DSGE.detexify("Φ") == "Phi"
    @test DSGE.detexify("Χ") == "Chi"
    @test DSGE.detexify("Ψ") == "Psi"
    @test DSGE.detexify("Ω") == "Omega"

    @test DSGE.detexify("′") == "'"
end

@testset "Check speye" begin
    @test DSGE.speye(2) == SparseMatrixCSC{Float64}(I, 2, 2)
    @test DSGE.speye(Complex, 2) == SparseMatrixCSC{Complex}(I, 2, 2)
end

@info "The following warning is expected."
@testset "Check macro for approximate equality of matrices" begin
    @test !(DSGE.test_matrix_eq2([0. 0.; 1. 2.] + [0. 0.; 1e-9 0.], [0. 0.; 1. 2.], "a", "a", 1e-10, 1e-10))
    @test DSGE.test_matrix_eq2([0. 0.; 1. 2.] + [0. 0.; 1e-9 0.], [0. 0.; 1. 2.], "a", "a", 1e-10, 2e-7)
end

@testset "Check overloading of operators for Complex numbers" begin
    @test Complex(1.,1.) > Complex(.99999, 1e10)
    @test 1. > Complex(.99999, 1e10)
    @test Complex(.99999, 1e10) < 1.
    @test min(Complex(1., 1.), Complex(.99999, 1e10)) == 0.99999
    @test min(1., Complex(.99999, 1e10)) == 0.99999
    @test min(Complex(1., 1.), 0.99999) == 0.99999
    @test max(Complex(1., 1.), Complex(.99999, 1e10)) == 1.
    @test max(1., Complex(.99999, 1e10)) == 1.
    @test max(Complex(1., 1.), 0.99999) == 1.
end

nothing
