using DSGE, JLD2, Dates, ModelConstructors, DataFrames
using Test, FileIO

path = dirname(@__FILE__)

# Initialize models
function make_test_model(year::Int)
    m = Model990()
    vint = Dates.format(Dates.Date(year, 4, 10), "yymmdd")
    m <= Setting(:data_vintage, vint)
    m <= Setting(:cond_vintage, vint)
    m <= Setting(:date_forecast_start, DSGE.quartertodate("$year-Q1"))
    m <= Setting(:date_conditional_end, DSGE.quartertodate("$year-Q1"))
    m <= Setting(:forecast_horizons, 16)
    return m
end

m_new = make_test_model(2016)
m_old = make_test_model(2014)

# Read in data and parameters
@load "$path/../reference/decompose_forecast_args.jld2" df_new df_old params_new params_old

# Read in expected outputs
@load "$path/../reference/decompose_forecast_out.jld2" exp_decomps

cond_types = [(:none, :none), (:none, :full), (:full, :none), (:full, :full)]
for (cond_new, cond_old) in cond_types
    # Test inner method
    decomps = decompose_forecast(m_new, m_old,
                                 cond_new == :full ? df_new : df_new[1:end-1, :],
                                 cond_old == :full ? df_old : df_old[1:end-1, :],
                                 params_new, params_old, cond_new, cond_old,
                                 [:obs, :pseudo]; check = true)

    for k in keys(decomps)
        @test decomps[k] ≈ exp_decomps[(cond_new, cond_old)][k]
    end
end

# Test outer method
m_new = AnSchorfheide(testing = true)
m_new <= Setting(:saveroot, tempdir())
m_new <= Setting(:date_forecast_start, quartertodate("2015-Q4"))
m_new <= Setting(:date_conditional_end, quartertodate("2015-Q4"))
m_new <= Setting(:use_population_forecast, true)
m_new <= Setting(:forecast_horizons, 12)

estroot = normpath(joinpath(dirname(@__FILE__), "..", "reference"))
overrides = forecast_input_file_overrides(m_new)
overrides[:mode] = joinpath(estroot, "optimize.h5")
overrides[:full] = joinpath(estroot, "mhsave_test.h5")
m_new <= Setting(:forecast_block_size, 20)

m_old = deepcopy(m_new)
m_old <= Setting(:date_forecast_start, quartertodate("2014-Q4"))
m_old <= Setting(:date_conditional_end, quartertodate("2014-Q4"))

df_new = load_data(m_new)
df_old = df_new[1:end-4, :]

decompose_forecast(m_new, m_old, df_new, df_old, :mode, :none, :none, [:obs, :pseudo];
                         verbose = :none)
decomposition_means(m_new, m_old, :mode, :none, :none, [:obs, :pseudo]; verbose = :none)
decompose_forecast(m_new, m_old, df_new, df_old, :full, :none, :none, [:obs, :pseudo], verbose = :none)
decomposition_means(m_new, m_old, :full, :none, :none, [:obs, :pseudo], verbose = :none)

## Regime switching
custom_settings = Dict{Symbol, Setting}(
    :data_vintage             => Setting(:data_vintage, "160812"),
    :cond_vintage             => Setting(:cond_vintage, "160812"),
    :cond_id                  => Setting(:cond_id, 0),
    :use_population_forecast  => Setting(:use_population_forecast, true),
    :date_presample_start     => Setting(:date_presample_start, Date(1959, 9, 30)),
    :date_forecast_start      => Setting(:date_forecast_start, DSGE.quartertodate("2016-Q3")),
    :date_conditional_end     => Setting(:date_conditional_end, DSGE.quartertodate("2016-Q3")),
    :forecast_horizons        => Setting(:forecast_horizons, 16),
    :n_mon_anticipated_shocks => Setting(:n_mon_anticipated_shocks, 6))
m    = Model1002("ss10", testing = true, custom_settings = custom_settings)  # baseline model
m_rs = Model1002("ss51v", testing = true, custom_settings = custom_settings) # pseudo regime switching (identical values for standard deviations)
m_rs <= Setting(:rate_expectations_source, :ois)
m_rs.settings[:regime_switching] = Setting(:regime_switching, true)
m_rs.settings[:n_regimes] = Setting(:n_regimes, 3)
m_rs.settings[:regime_switching_ndraws] = Setting(:regime_switching_ndraws, 4)
m_rs.test_settings[:regime_switching] = Setting(:regime_switching, true)
m_rs.test_settings[:n_regimes] = Setting(:n_regimes, 3)
m_rs.test_settings[:regime_switching_ndraws] = Setting(:regime_switching_ndraws, 4)
m_rs.settings[:regime_dates] = Setting(:regime_dates,
                                       Dict{Int, Date}(1 => date_presample_start(m), 2 => Date(2010, 3, 31), 3 => Date(2012, 9, 30)))
m_rs.test_settings[:regime_dates] = Setting(:regime_dates,
                                            Dict{Int, Date}(1 => date_presample_start(m), 2 => Date(2010, 3, 31), 3 => Date(2012, 9, 30)))
setup_regime_switching_inds!(m_rs)
df = load("$path/../reference/regime_switch_data.jld2", "none")

for i in 1:3
    adj = (i == 1) ? 1. : .95
    ModelConstructors.set_regime_val!(m_rs[:α], i, adj * m[:α].value; override_bounds = true)
    ModelConstructors.set_regime_val!(m_rs[:σ_g], i, adj * m[:σ_g].value; override_bounds = true)
    ModelConstructors.set_regime_val!(m_rs[:σ_b], i, adj * m[:σ_b].value; override_bounds = true)
    ModelConstructors.set_regime_val!(m_rs[:σ_μ], i, adj * m[:σ_μ].value; override_bounds = true)
    ModelConstructors.set_regime_val!(m_rs[:σ_ztil], i, adj * m[:σ_ztil].value; override_bounds = true)
    ModelConstructors.set_regime_val!(m_rs[:σ_λ_f], i, adj * m[:σ_λ_f].value; override_bounds = true)
    ModelConstructors.set_regime_val!(m_rs[:σ_λ_w], i, adj * m[:σ_λ_w].value; override_bounds = true)
    ModelConstructors.set_regime_val!(m_rs[:σ_r_m], i, adj * m[:σ_r_m].value; override_bounds = true)
    ModelConstructors.set_regime_val!(m_rs[:σ_σ_ω], i, adj * m[:σ_σ_ω].value; override_bounds = true)
    ModelConstructors.set_regime_val!(m_rs[:σ_μ_e], i, adj * m[:σ_μ_e].value; override_bounds = true)
    ModelConstructors.set_regime_val!(m_rs[:σ_γ], i, adj * m[:σ_γ].value; override_bounds = true)
    ModelConstructors.set_regime_val!(m_rs[:σ_π_star], i, adj * m[:σ_π_star].value; override_bounds = true)
    ModelConstructors.set_regime_val!(m_rs[:σ_lr], i, adj * m[:σ_lr].value; override_bounds = true)
    ModelConstructors.set_regime_val!(m_rs[:σ_z_p], i, adj * m[:σ_z_p].value; override_bounds = true)
    ModelConstructors.set_regime_val!(m_rs[:σ_tfp], i, adj * m[:σ_tfp].value; override_bounds = true)
    ModelConstructors.set_regime_val!(m_rs[:σ_gdpdef], i, adj * m[:σ_gdpdef].value; override_bounds = true)
    ModelConstructors.set_regime_val!(m_rs[:σ_corepce], i, adj * m[:σ_corepce].value; override_bounds = true)
    ModelConstructors.set_regime_val!(m_rs[:σ_gdp], i, adj * m[:σ_gdp].value; override_bounds = true)
    ModelConstructors.set_regime_val!(m_rs[:σ_gdi], i, adj * m[:σ_gdi].value; override_bounds = true)

    for j = 1:DSGE.n_mon_anticipated_shocks(m_rs)
        ModelConstructors.set_regime_val!(m_rs[Symbol("σ_r_m$(j)")], i, adj * m[Symbol("σ_r_m$(j)")]; override_bounds = true)
    end
end

m_rs_old = deepcopy(m_rs)
m_rs_old <= Setting(:date_forecast_start, quartertodate("2015-Q3"))
m_rs_old <= Setting(:date_conditional_end, quartertodate("2015-Q3"))
m_rs_old.settings[:forecast_horizons] = Setting(:forecast_horizons, 12)
m_rs_old.test_settings[:forecast_horizons] = Setting(:forecast_horizons, 12)
df_old = df[1:end - 4, :]

decompose_forecast(m_rs, m_rs_old, df, df_old, map(x -> x.value, m_rs.parameters),
                   map(x -> x.value, m_rs_old.parameters), :none, :none, [:obs, :pseudo])
@test_broken decomposition_means(m_rs, m_rs_old, :mode, :none, :none, [:obs, :pseudo]; verbose = :none) # TODO: NEED METHOD THAT LOADS THE DRAWS FROM A SAVED FILE B/C THAT ONE WRITES TO STUFF TO A JLD2 FILE

nothing
