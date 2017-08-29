using DSGE, DataFrames, JLD

path = dirname(@__FILE__())

# Set up arguments
m = AnSchorfheide(testing = true)
m <= Setting(:date_forecast_start, quartertodate("2015-Q4"))

df, system, kal = jldopen("$path/../reference/forecast_args.jld","r") do file
    read(file, "df"), read(file, "system"), read(file, "kalman")
end

# Read expected output
exp_states, exp_shocks, exp_pseudo = jldopen("$path/../reference/smooth_out.jld", "r") do file
    read(file, "exp_states"),
    read(file, "exp_shocks"),
    read(file, "exp_pseudo")
end

# Smooth without drawing states
states = Dict{Symbol, Matrix{Float64}}()
shocks = Dict{Symbol, Matrix{Float64}}()
pseudo = Dict{Symbol, Matrix{Float64}}()

for smoother in [:hamilton, :koopman, :carter_kohn, :durbin_koopman]
    m <= Setting(:forecast_smoother, smoother)

    states[smoother], shocks[smoother], pseudo[smoother] =
        smooth(m, df, system, kal; draw_states = false)

    @test_matrix_approx_eq exp_states states[smoother]
    @test_matrix_approx_eq exp_shocks shocks[smoother]
    @test_matrix_approx_eq exp_pseudo pseudo[smoother]
end

# Smooth, drawing states
for smoother in [:carter_kohn, :durbin_koopman]
    m <= Setting(:forecast_smoother, smoother)
    smooth(m, df, system, kal; draw_states = true)
end


nothing
