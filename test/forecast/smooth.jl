path = dirname(@__FILE__())

# Set up arguments
m = AnSchorfheide(testing = true)
m <= Setting(:date_forecast_start, quartertodate("2015-Q4"))

forecast_args = load("$path/../reference/forecast_args.jld2")
df = forecast_args["df"]
system = forecast_args["system"]

# Read expected output
smooth_out = load("$path/../reference/smooth_out.jld2")
exp_states = smooth_out["exp_states"]
exp_shocks = smooth_out["exp_shocks"]
exp_pseudo = smooth_out["exp_pseudo"]

# Smooth without drawing states
states = Dict{Symbol, Matrix{Float64}}()
shocks = Dict{Symbol, Matrix{Float64}}()
pseudo = Dict{Symbol, Matrix{Float64}}()

@testset "Test smoother without drawing states" begin
    for smoother in [:hamilton, :koopman, :carter_kohn, :durbin_koopman]
        m <= Setting(:forecast_smoother, smoother)

        states[smoother], shocks[smoother], pseudo[smoother] =
            smooth(m, df, system; draw_states = false)

        @test @test_matrix_approx_eq exp_states states[smoother]
        @test @test_matrix_approx_eq exp_shocks shocks[smoother]
        @test @test_matrix_approx_eq exp_pseudo pseudo[smoother]
    end
end

# Smooth, drawing states
for smoother in [:carter_kohn, :durbin_koopman]
    m <= Setting(:forecast_smoother, smoother)
    smooth(m, df, system; draw_states = true)
end


nothing
