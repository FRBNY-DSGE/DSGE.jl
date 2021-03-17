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

# Check trends for "fake" regime-switching system
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
m.settings[:regime_switching] = Setting(:regime_switching, true)
m.test_settings[:regime_switching] = Setting(:regime_switching, true)
setup_regime_switching_inds!(m)
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

## Shock decompositions with time-varying CCC
# Set up
m = Model1002("ss10"; custom_settings = Dict{Symbol, Setting}(:add_altpolicy_pgap => Setting(:add_altpolicy_pgap, true),
                                                              :add_altpolicy_ygap => Setting(:add_altpolicy_ygap, true)))
m <= Setting(:regime_switching, true)
m <= Setting(:regime_dates, Dict{Int, Date}(1 => date_presample_start(m),
                                            2 => Date(2020, 6, 30),
                                            3 => Date(2020, 9, 30),
                                            4 => Date(2020, 12, 31),
                                            5 => Date(2021, 3, 31),
                                            6 => Date(2021, 6, 30)))
m <= Setting(:date_forecast_start, Date(2020, 6, 30))
m <= Setting(:date_conditional_end, Date(2020, 6, 30))
m <= Setting(:tvis_information_set, [1:1, 2:2, 3:6, 4:6, 5:6, 6:6])
m <= Setting(:replace_eqcond, true)
m <= Setting(:gensys2, true)
zlb_rule_eqcond = DSGE.EqcondEntry(DSGE.zlb_rule(), [1., 0.])
m <= Setting(:regime_eqcond_info, Dict(3 => deepcopy(zlb_rule_eqcond),
                                             4 => deepcopy(zlb_rule_eqcond),
                                             5 => deepcopy(zlb_rule_eqcond),
                                             6 => DSGE.EqcondEntry(DSGE.flexible_ait(), [1., 0.])))
m <= Setting(:temporary_altpolicy_names, [:zero_rate])
setup_regime_switching_inds!(m; cond_type = :full)
df = load(joinpath(path, "..", "reference", "regime_switch_data.jld2"), "regime_switch_df_full")
sys = compute_system(m; tvis = true)
_, histshocks, _, init_states = smooth(m, df, sys; cond_type = :full)
output = DSGE.forecast_one_draw(m, :mode, :full, [:forecastobs, :histpseudo, :forecastpseudo,
                                                  :histstates, :forecaststates],
                                [x.value for x in m.parameters], df, regime_switching = true,
                                n_regimes = get_setting(m, :n_regimes))
shockstates, shockobs, shockpseudo = shock_decompositions(m, sys, histshocks,
                                                          date_mainsample_start(m), date_conditional_end(m), :full)
dettrendstates, dettrendobs, dettrendpseudo = deterministic_trends(m, sys, init_states,
                                                                   date_mainsample_start(m),
                                                                   date_conditional_end(m), :full)
trendstates, trendobs, trendpseudo = trends(m, sys, date_mainsample_start(m), date_conditional_end(m), :full)

states = hcat(output[:histstates], output[:forecaststates])[:, index_shockdec_start(m):end]
obs = hcat(df_to_matrix(m, df; cond_type = :full, include_presample = false)[:, 1:end - 1],
           output[:forecastobs])[:, index_shockdec_start(m):end]
pseudo = hcat(output[:histpseudo], output[:forecastpseudo])[:, index_shockdec_start(m):end]

impl_states = dropdims(sum(shockstates, dims = 3), dims = 3) + dettrendstates + trendstates
impl_obs = dropdims(sum(shockobs, dims = 3), dims = 3) + dettrendobs + trendobs
impl_pseudo = dropdims(sum(shockpseudo, dims = 3), dims = 3) + dettrendpseudo + trendpseudo

@test states ≈ impl_states
@test obs[vcat(1:9, 13), :] ≈ impl_obs[vcat(1:9, 13), :]
@test pseudo ≈ impl_pseudo

df_states = DSGE.prepare_means_table_trend_nostates(m, :full, :state,
                                                    date_mainsample_start(m), date_conditional_end(m), apply_altpolicy = true,
                                                    annualize = false)
df_obs    = DSGE.prepare_means_table_trend_nostates(m, :full, :obs,
                                                    date_mainsample_start(m), date_conditional_end(m), apply_altpolicy = true,
                                                    annualize = false)
df_pseudo = DSGE.prepare_means_table_trend_nostates(m, :full, :pseudo,
                                                    date_mainsample_start(m), date_conditional_end(m), apply_altpolicy = true,
                                                    annualize = false)

@test all(Matrix(df_states[:, 2:end]) .≈ 0.)
@test all((df_to_matrix(m, df_obs) .- sys[n_regimes(sys), :DD]) .≈ 0.)
pseudo_obs = vcat(1:13, 19:n_pseudo_observables(m)) .+ 1 # Remove forward-looking pseudo-obs
@test all((Matrix(df_pseudo[:, pseudo_obs])' .- sys[n_regimes(sys), :DD_pseudo][pseudo_obs .- 1]) .≈ 0.)

nothing
