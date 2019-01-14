using DSGE
using Base.Test, BenchmarkTools
using JLD

# What do you want to do?
check_steady_state = true
check_jacobian = true
check_solution = true

path = dirname(@__FILE__)

m = RealBond()

# Steady-state computation
if check_steady_state
    steadystate!(m)
    # @btime steadystate!(m)

    file = jldopen("$path/reference/steady_state.jld", "r")
    saved_ell  = read(file, "ell")
    saved_c    = read(file, "c")
    saved_η    = read(file, "eta")
    saved_μ    = read(file, "mu")
    saved_β    = read(file, "beta")
    saved_χss  = read(file, "chi_ss")
    close(file)

    @testset "Check steady state outputs" begin
        @test saved_ell ≈ m[:lstar].value
        @test saved_c   ≈ m[:cstar].value
        @test saved_η   ≈ m[:ηstar].value
        @test saved_μ   ≈ m[:μstar].value
        @test saved_β   ≈ m[:βstar].value
        @test saved_χss ≈ m[:χstar].value
    end
end

# Jacobian computation
if check_jacobian
    m.testing = true        # So that it will test against the unnormalized Jacobian
    JJ = DSGE.jacobian(m)
    # @btime JJ = DSGE.jacobian(m)

    file = jldopen("$path/reference/jacobian.jld", "r")
    saved_JJ  = read(file, "JJ")

    MUP = read(file, "MUP")
    ZP  = read(file, "ZP")
    MONP = read(file, "MONP")
    ELLP = read(file, "ELLP")
    RRP = read(file, "RRP")
    IIP = read(file, "IIP")
    WWP = read(file, "WWP")
    PIP = read(file, "PIP")
    TTP = read(file, "TTP")
    MU  = read(file, "MU")
    Z   = read(file, "Z")
    MON = read(file, "MON")
    ELL = read(file, "ELL")
    RR  = read(file, "RR")
    II  = read(file, "II")
    WW  = read(file, "WW")
    PI  = read(file, "PI")
    TT  = read(file, "TT")

    F1 = read(file, "F1")
    F2 = read(file, "F2")
    F3 = read(file, "F3")
    F4 = read(file, "F4")
    F5 = read(file, "F5")
    F6 = read(file, "F6")
    F7 = read(file, "F7")
    F8 = read(file, "F8")
    F9 = read(file, "F9")
    close(file)

    @testset "Check jacobian outputs" begin
        @testset "Euler Equation" begin
            @test saved_JJ[F1, ELLP] ≈ JJ[F1, ELLP]
            @test saved_JJ[F1, RRP]  ≈ JJ[F1, RRP]
            @test saved_JJ[F1, WWP]  ≈ JJ[F1, WWP]
            @test saved_JJ[F1, TTP]  ≈ JJ[F1, TTP]
            @test saved_JJ[F1, ELL]  ≈ JJ[F1, ELL]
            @test saved_JJ[F1, RR]   ≈ JJ[F1, RR]
            @test saved_JJ[F1, WW]   ≈ JJ[F1, WW]
            @test saved_JJ[F1, TT]   ≈ JJ[F1, TT]
        end

        @testset "Kolmogorov Forward Equation" begin
            @test saved_JJ[F2, MUP] ≈ JJ[F2, MUP]
            @test saved_JJ[F2, MU]  ≈ JJ[F2, MU]
            @test saved_JJ[F2, ELL] ≈ JJ[F2, ELL]
            @test saved_JJ[F2, WW]  ≈ JJ[F2, WW]
            @test saved_JJ[F2, RR]  ≈ JJ[F2, RR]
            @test saved_JJ[F2, TT]  ≈ JJ[F2, TT]
        end

        @testset "Market clearing condition" begin
            @test saved_JJ[F3, MU]  ≈ JJ[F3, MU]
            @test saved_JJ[F3, Z]   ≈ JJ[F3, Z]
            @test saved_JJ[F3, ELL] ≈ JJ[F3, ELL]
            @test saved_JJ[F3, WW]  ≈ JJ[F3, WW]
            @test saved_JJ[F3, RR]  ≈ JJ[F3, RR]
            @test saved_JJ[F3, TT]  ≈ JJ[F3, TT]
        end

        @testset "Technology process" begin
            @test saved_JJ[F4, ZP] ≈ JJ[F4, ZP]
            @test saved_JJ[F4, Z] ≈ JJ[F4, Z]
        end

        @testset "Phillips equation" begin
            @test saved_JJ[F5, PIP] ≈ JJ[F5, PIP]
            @test saved_JJ[F5, Z]   ≈ JJ[F5, Z]
            @test saved_JJ[F5, WW]  ≈ JJ[F5, WW]
            @test saved_JJ[F5, PI]  ≈ JJ[F5, PI]
        end

        @testset "Taylor equation" begin
            @test saved_JJ[F6, II]  ≈ JJ[F6, II]
            @test saved_JJ[F6, PI]  ≈ JJ[F6, PI]
            @test saved_JJ[F6, MON] ≈ JJ[F6, MON]
        end

        @testset "Fisher equation" begin
            @test saved_JJ[F7, II]  ≈ JJ[F7, II]
            @test saved_JJ[F7, PIP] ≈ JJ[F7, PIP]
            @test saved_JJ[F7, RR]  ≈ JJ[F7, RR]
        end

        @testset "Transfers equation" begin
            @test saved_JJ[F8, TT] ≈ JJ[F8, TT]
            @test saved_JJ[F8, Z]  ≈ JJ[F8, Z]
            @test saved_JJ[F8, WW] ≈ JJ[F8, WW]
        end
    end
end

# Solve
if check_solution
    m.testing = false      # So the Jacobian will be normalized within the klein solution
    gx, hx = klein(m)
    # @btime klein(m)

    file = jldopen("$path/reference/solve.jld", "r")
    saved_gx   = read(file, "gx")
    saved_hx   = read(file, "hx")
    close(file)

    @testset "Check solve outputs" begin
        @test saved_gx  ≈ gx
        @test saved_hx  ≈ hx
    end
end
