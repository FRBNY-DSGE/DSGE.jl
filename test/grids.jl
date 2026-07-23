using DSGE, Test, BenchmarkTools

# Minimal concrete model so get_grid can be exercised without constructing a
# full (expensive) heterogeneous-agent model — get_grid only touches `m.grids`.
struct GridTestModel{T} <: DSGE.AbstractDSGEModel{T}
    grids::Dict{Symbol, Grid}
end

@testset "grids.jl" begin

    @testset "Grid type and constructors" begin
        points  = collect(range(0.0, stop = 1.0, length = 5))
        weights = fill(0.2, 5)

        g = Grid(points, weights, 1.0)
        @test g.points  == points
        @test g.weights == weights
        @test g.scale   == 1.0

        # Inner constructor rejects weights that don't sum to `scale`
        @test_throws ErrorException Grid(points, fill(0.5, 5), 1.0)

        # Quadrature-based constructor (uniform): weights sum to scale (default 1)
        gu = Grid(DSGE.uniform_quadrature, 0.0, 1.0, 5)
        @test length(gu.points) == 5
        @test sum(gu.weights) ≈ 1.0
        @test gu.points[1] ≈ 0.0 && gu.points[end] ≈ 1.0
    end

    @testset "uniform_quadrature" begin
        grid, wts = DSGE.uniform_quadrature(0.0, 1.0, 10)
        @test length(grid) == 10
        @test length(wts)  == 10
        @test grid ≈ collect(range(0.0, stop = 1.0, length = 10))
        @test sum(wts) ≈ 1.0
        @test all(wts .≈ 1 / 10)

        # Curried form bakes in the scale
        q = DSGE.uniform_quadrature(2.0)
        grid2, wts2 = q(0.0, 1.0, 10)
        @test grid2 ≈ grid
        @test sum(wts2) ≈ 2.0
    end

    # curtis_clenshaw_quadrature depends on `chebpts`, which lives in chebyshev.jl.
    # That file is NOT included by DSGE.jl, so chebpts is undefined and this rule is
    # currently dead code. Guard on availability so the suite passes now and these
    # tests light up automatically if chebyshev.jl is ever wired into the package.
    @testset "curtis_clenshaw_quadrature" begin
        if isdefined(DSGE, :chebpts)
            g, w = DSGE.curtis_clenshaw_quadrature(-1.0, 1.0, 9)
            @test length(vec(g)) == 9
            @test length(vec(w)) == 9
            @test all(-1.0 - 1e-10 .<= vec(g) .<= 1.0 + 1e-10)
            @test sum(vec(w)) ≈ 2.0    # ∫_{-1}^{1} 1 dx = interval length

            # On a shifted interval the weights still integrate to the interval length
            g2, w2 = DSGE.curtis_clenshaw_quadrature(0.0, 4.0, 9)
            @test sum(vec(w2)) ≈ 4.0
            @test all(0.0 - 1e-10 .<= vec(g2) .<= 4.0 + 1e-10)

            # Grid constructor via the Clenshaw-Curtis rule (scale = interval length)
            gc = Grid(DSGE.curtis_clenshaw_quadrature, -1.0, 1.0, 7; scale = 2.0)
            @test length(gc.points) == 7
            @test sum(gc.weights) ≈ 2.0

            # Curried form bakes in the Chebyshev kind
            q = DSGE.curtis_clenshaw_quadrature(2)
            g3, w3 = q(-1.0, 1.0, 9)
            @test vec(g3) ≈ vec(g)
            @test vec(w3) ≈ vec(w)
        else
            @info "Skipping curtis_clenshaw_quadrature tests: chebpts (chebyshev.jl) is not loaded"
            @test_skip DSGE.curtis_clenshaw_quadrature(-1.0, 1.0, 9)
        end
    end

    @testset "tauchen86 (AR(1))" begin
        μ, ρ, σ, n, λ = 0.1, 0.8, 0.2, 7, 3.0
        xgrid, xprob, xscale = DSGE.tauchen86(μ, ρ, σ, n, λ)

        @test length(xgrid) == n
        @test size(xprob)   == (n, n)

        # Transition matrix rows are valid probability distributions
        @test vec(sum(xprob, dims = 2)) ≈ ones(n)
        @test all(xprob .>= 0.0)

        # Grid is evenly spaced by xscale and centered at the unconditional mean
        center = μ / (1 - ρ)
        @test xscale ≈ 2λ * σ / ((1 - ρ) * (n - 1))
        @test diff(xgrid) ≈ fill(xscale, n - 1)
        @test (xgrid[1] + xgrid[end]) / 2 ≈ center
    end

    @testset "tauchen86 (i.i.d.)" begin
        μ, σ, n, λ = 0.1, 0.2, 7, 3.0
        xg, p, xs    = DSGE.tauchen86(μ, σ, n, λ)
        xgA, xpA, xsA = DSGE.tauchen86(μ, 0.0, σ, n, λ)

        @test length(xg) == n
        @test length(p)  == n
        # The 4-arg method delegates to the 5-arg method with ρ = 0
        @test xg ≈ xgA
        @test xs ≈ xsA
        @test p  ≈ vec(xpA[:, 1])
    end

    @testset "get_grid" begin
        g  = Grid(collect(range(0.0, stop = 1.0, length = 5)), fill(0.2, 5), 1.0)
        gm = GridTestModel{Float64}(Dict(:test_grid => g))
        @test DSGE.get_grid(gm, :test_grid) === g
        @test isa(DSGE.get_grid(gm, :test_grid), Grid)
        @test_throws KeyError DSGE.get_grid(gm, :missing_grid)
    end

    @testset "quadrature_sum" begin
        g = Grid(collect(range(0.0, stop = 1.0, length = 5)), fill(0.2, 5), 1.0)
        x = [1.0, 2.0, 3.0, 4.0, 5.0]
        @test DSGE.quadrature_sum(x, g) ≈ sum(g.weights .* x .* g.points)

        # NOTE: the (grid, x) argument order calls `sum(x, grid)` internally, which
        # dispatches to `sum(f, itr)` with a non-callable `x` and errors. Marked
        # broken; if the overload is fixed this will flag as an unexpected pass.
        @test_broken DSGE.quadrature_sum(g, x) ≈ DSGE.quadrature_sum(x, g)
    end
end

####################
# Benchmarks
####################
@info "Benchmarking grids.jl"
let
    pts = collect(range(0.0, stop = 1.0, length = 100))
    wts = fill(0.01, 100)
    g   = Grid(pts, wts, 1.0)
    x   = collect(1.0:100.0)

    @btime Grid($pts, $wts, 1.0)
    @btime DSGE.uniform_quadrature(0.0, 1.0, 100)
    isdefined(DSGE, :chebpts) && @btime DSGE.curtis_clenshaw_quadrature(-1.0, 1.0, 100)
    @btime DSGE.tauchen86(0.1, 0.8, 0.2, 50, 3.0)
    @btime DSGE.tauchen86(0.1, 0.2, 50, 3.0)
    @btime DSGE.quadrature_sum($x, $g)
end

nothing
