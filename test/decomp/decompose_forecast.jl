using DSGE, JLD2
using Test

path = dirname(@__FILE__)

# Initialize models
function make_test_model(year::Int)
    m = Model990()
    vint = Dates.format(Date(year, 4, 10), "yymmdd")
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
    @time decomps = decompose_forecast(m_new, m_old,
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

m_old = deepcopy(m_new)
m_old <= Setting(:date_forecast_start, quartertodate("2014-Q4"))
m_old <= Setting(:date_conditional_end, quartertodate("2014-Q4"))

df_new = load_data(m_new)
df_old = df_new[1:end-4, :]

@time decompose_forecast(m_new, m_old, df_new, df_old, :mode, :none, :none, [:obs, :pseudo];
                         verbose = :none)
@time decomposition_means(m_new, m_old, :mode, :none, :none, [:obs, :pseudo], verbose = :none)


nothing
