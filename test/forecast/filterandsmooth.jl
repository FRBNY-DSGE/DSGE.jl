using DSGE, DataFrames, JLD
include("../util.jl")

path = dirname(@__FILE__())

# Set up arguments
m = AnSchorfheide(testing = true)
m <= Setting(:date_forecast_start, quartertodate("2015-Q4"))

df, system, z0, P0 = jldopen("$path/../reference/forecast_args.jld","r") do file
    read(file, "df"), read(file, "system"), read(file, "z0"), read(file, "vz0")
end

# Read in expected output
exp_states, exp_shocks, exp_pseudo =
    jldopen("$path/../reference/smooth_out.jld", "r") do file
        read(file, "exp_states"),
        read(file, "exp_shocks"),
        read(file, "exp_pseudo")
    end

for smoother in [:durbin_koopman, :koopman]
    m <= Setting(:forecast_smoother, smoother)

    # Temporary (during renaming)
    if smoother == :koopman
        smoother = :kalman
    end

    # Without providing z0 and P0
    states, shocks, pseudo = filterandsmooth(m, df, system)

    @test_matrix_approx_eq exp_states[(smoother, :no_z0)] states
    @test_matrix_approx_eq exp_shocks[(smoother, :no_z0)] shocks
    @test_matrix_approx_eq exp_pseudo[(smoother, :no_z0)] pseudo

    # Providing z0 and P0
    states, shocks, pseudo = filterandsmooth(m, df, system, z0, P0)

    @test_matrix_approx_eq exp_states[(smoother, :z0)] states
    @test_matrix_approx_eq exp_shocks[(smoother, :z0)] shocks
    @test_matrix_approx_eq exp_pseudo[(smoother, :z0)] pseudo
end

nothing