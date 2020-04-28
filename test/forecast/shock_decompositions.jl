using DSGE, FileIO, JLD2, ModelConstructors, Test, Random, Dates
path = dirname(@__FILE__)

# Set up arguments
m = AnSchorfheide(testing = true)
m <= Setting(:date_forecast_start, quartertodate("2015-Q4"))
m <= Setting(:forecast_horizons, 1)

system, histshocks = JLD2.jldopen("$path/../reference/forecast_args.jld2","r") do file
    read(file, "system"), read(file, "histshocks")
end

# Read expected output
exp_states, exp_obs, exp_pseudo =
    JLD2.jldopen("$path/../reference/shock_decompositions_out.jld2", "r") do file
        read(file, "exp_states"), read(file, "exp_obs"), read(file, "exp_pseudo")
    end

# With shockdec_startdate not null
states, obs, pseudo = shock_decompositions(m, system, histshocks)

@testset "Test shockdec with non-null startdate" begin
    @test @test_matrix_approx_eq exp_states[:startdate] states
    @test @test_matrix_approx_eq exp_obs[:startdate]    obs
    @test @test_matrix_approx_eq exp_pseudo[:startdate] pseudo
end

# With shockdec_startdate null
#m <= Setting(:shockdec_startdate, Nullable{Date}())
m <= Setting(:shockdec_startdate, nothing)
states, obs, pseudo = shock_decompositions(m, system, histshocks)

@testset "Test shockdec with null startdate" begin
    @test @test_matrix_approx_eq exp_states[:no_startdate] states
    @test @test_matrix_approx_eq exp_obs[:no_startdate]    obs
    @test @test_matrix_approx_eq exp_pseudo[:no_startdate] pseudo
end

@testset "Deterministic trends" begin
    dettrend = DSGE.deterministic_trends(m, system, zeros(8))
    @test @test_matrix_approx_eq dettrend[1] zeros(size(dettrend[1]))
    @test @test_matrix_approx_eq dettrend[2] zeros(size(dettrend[2]))
    @test @test_matrix_approx_eq dettrend[3] zeros(size(dettrend[3]))
end

@testset "Trends" begin
    out = DSGE.trends(system)
    @test @test_matrix_approx_eq out[1] system[:CCC]
    @test @test_matrix_approx_eq out[2] system[:ZZ] * system[:CCC] + system[:DD]
    @test @test_matrix_approx_eq out[3] system[:ZZ_pseudo] * system[:CCC] + system[:DD_pseudo]
end

# Check trends for regime-switching system
reg_sys = RegimeSwitchingSystem(System{Float64}[system, system])
@testset "Regime-switching trends" begin
    out = DSGE.trends(reg_sys)
    for i in 1:2
        @test @test_matrix_approx_eq out[1][:, i] reg_sys[i, :CCC]
        @test @test_matrix_approx_eq out[2][:, i] reg_sys[i, :ZZ] * reg_sys[i, :CCC] + reg_sys[i, :DD]
        @test @test_matrix_approx_eq out[3][:, i] reg_sys[i, :ZZ_pseudo] * reg_sys[i, :CCC] + reg_sys[i, :DD_pseudo]
    end
end

m.settings[:regime_dates] = Setting(:regime_dates, Dict{Int, Date}(1 => date_presample_start(m), 2 => Date(2010, 3, 31)))
m.test_settings[:regime_dates] = Setting(:regime_dates, Dict{Int, Date}(1 => date_presample_start(m), 2 => Date(2010, 3, 31)))
@testset "Regime-switching deterministic trends" begin
    dettrend = DSGE.deterministic_trends(m, reg_sys, zeros(8))
    @test @test_matrix_approx_eq dettrend[1] zeros(size(dettrend[1]))
    @test @test_matrix_approx_eq dettrend[2] zeros(size(dettrend[2]))
    @test @test_matrix_approx_eq dettrend[3] zeros(size(dettrend[3]))
end


# Check shock decompositions for regime-switching system
out_shockdec1 = shock_decompositions(m, reg_sys, histshocks,
                                    date_presample_start(m), date_forecast_start(m)) # check we only return the shockdecs for requested period
out_shockdec2 = shock_decompositions(m, reg_sys, histshocks[:, 1:end - 1]) # check default

@testset "Shock decompositions with regime-switching" begin
    @test @test_matrix_approx_eq out_shockdec1[1] states
    @test @test_matrix_approx_eq out_shockdec1[2] obs
    @test @test_matrix_approx_eq out_shockdec1[3] pseudo
    @test @test_matrix_approx_eq out_shockdec2[1] states
    @test @test_matrix_approx_eq out_shockdec2[2] obs
    @test @test_matrix_approx_eq out_shockdec2[3] pseudo
end

nothing
