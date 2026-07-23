using DSGE, Test, BenchmarkTools
using Distributions: Normal, pdf, cdf

# Bring the (unexported) financial-frictions helpers into scope tersely
const FF = DSGE

# Central finite difference
cdiff(f, x; h = 1e-6) = (f(x + h) - f(x - h)) / (2h)

# A representative steady-state point (z = standardized default threshold, σ = idiosyncratic
# vol, spr = gross credit spread). These values yield sane nk, μ ∈ (0,1).
z0, σ0, spr0 = 0.3, 0.5, 1.01

@testset "financial_frictions.jl" begin

    @testset "Base functions and identities" begin
        # Closed forms
        @test FF.ω_fn(z0, σ0) ≈ exp(σ0 * z0 - σ0^2 / 2)
        @test FF.G_fn(z0, σ0) ≈ cdf(Normal(), z0 - σ0)
        @test FF.Γ_fn(z0, σ0) ≈ FF.ω_fn(z0, σ0) * (1 - cdf(Normal(), z0)) + cdf(Normal(), z0 - σ0)

        # Log-normal identity that underlies the derivative algebra: ω̄·φ(z) = φ(z−σ)
        @test FF.ω_fn(z0, σ0) * pdf(Normal(), z0) ≈ pdf(Normal(), z0 - σ0)
    end

    @testset "ω̄-derivatives vs. numerical (vary z at fixed σ, chain through ω̄)" begin
        dω̄_dz = cdiff(z -> FF.ω_fn(z, σ0), z0)

        # First derivatives: dX/dω̄ = (dX/dz) / (dω̄/dz)
        @test FF.dΓ_dω_fn(z0)      ≈ cdiff(z -> FF.Γ_fn(z, σ0), z0) / dω̄_dz  rtol = 1e-5
        @test FF.dG_dω_fn(z0, σ0)  ≈ cdiff(z -> FF.G_fn(z, σ0), z0) / dω̄_dz  rtol = 1e-5

        # Second derivatives: differentiate the analytic first derivative w.r.t. ω̄
        @test FF.d2Γ_dω2_fn(z0, σ0) ≈ cdiff(z -> FF.dΓ_dω_fn(z),     z0) / dω̄_dz  rtol = 1e-5
        @test FF.d2G_dω2_fn(z0, σ0) ≈ cdiff(z -> FF.dG_dω_fn(z, σ0), z0) / dω̄_dz  rtol = 1e-5
    end

    @testset "σ-derivatives vs. numerical (vary σ holding ω̄ fixed)" begin
        # Holding ω̄ fixed means z must be recomputed as σ moves: z(σ) = (ln ω̄ + σ²/2)/σ
        logω̄ = log(FF.ω_fn(z0, σ0))
        zfun(σ) = (logω̄ + σ^2 / 2) / σ
        @test zfun(σ0) ≈ z0   # sanity: the reparameterization recovers z0 at σ0

        @test FF.dG_dσ_fn(z0, σ0) ≈ cdiff(σ -> FF.G_fn(zfun(σ), σ), σ0)  rtol = 1e-5
        @test FF.dΓ_dσ_fn(z0, σ0) ≈ cdiff(σ -> FF.Γ_fn(zfun(σ), σ), σ0)  rtol = 1e-5

        # Cross-partials: ∂(dX/dω̄)/∂σ holding ω̄ fixed
        @test FF.d2G_dωdσ_fn(z0, σ0) ≈ cdiff(σ -> FF.dG_dω_fn(zfun(σ), σ), σ0)  rtol = 1e-5
        @test FF.d2Γ_dωdσ_fn(z0, σ0) ≈ cdiff(σ -> FF.dΓ_dω_fn(zfun(σ)),    σ0)  rtol = 1e-5
    end

    @testset "Spread elasticities and steady-state ratios" begin
        nk = FF.nk_fn(z0, σ0, spr0)
        μ  = FF.μ_fn(z0, σ0, spr0)

        # Net-worth-to-capital and monitoring-cost fractions must be sensible
        @test 0 < nk < 1
        @test 0 < μ  < 1

        # Elasticities are finite real numbers at the steady-state point
        @test isfinite(FF.ζ_zω_fn(z0, σ0, spr0))
        @test isfinite(FF.ζ_bω_fn(z0, σ0, spr0))
        @test isfinite(FF.ζ_spb_fn(z0, σ0, spr0))

        # ζ_spb_fn is the documented composition of ζ_bω/ζ_zω, nk
        zr = FF.ζ_bω_fn(z0, σ0, spr0) / FF.ζ_zω_fn(z0, σ0, spr0)
        @test FF.ζ_spb_fn(z0, σ0, spr0) ≈ -zr / (1 - zr) * nk / (1 - nk)
    end
end

####################
# Benchmarks
####################
run_benchmarks = false

if run_benchmarks
    @info "Benchmarking financial_frictions.jl"
    let
        @btime FF.ω_fn($z0, $σ0)
        @btime FF.G_fn($z0, $σ0)
        @btime FF.Γ_fn($z0, $σ0)
        @btime FF.μ_fn($z0, $σ0, $spr0)
        @btime FF.nk_fn($z0, $σ0, $spr0)
        @btime FF.ζ_zω_fn($z0, $σ0, $spr0)
        @btime FF.ζ_bω_fn($z0, $σ0, $spr0)
        @btime FF.ζ_spb_fn($z0, $σ0, $spr0)
    end
end

nothing
