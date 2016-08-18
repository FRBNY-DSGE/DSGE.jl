"""
```
default_settings!(m::AbstractModel)
```

Default Settings are constructed, initialized and added to `m.settings`.
"""
function default_settings!(m::AbstractModel)

    # I/O File locations
    saveroot = normpath(joinpath(dirname(@__FILE__), "..","save"))
    datapath = normpath(joinpath(dirname(@__FILE__), "..","save","input_data"))

    m <= Setting(:saveroot, saveroot, "Root of data directory structure")
    m <= Setting(:dataroot, datapath, "Input data directory path")

    # Data settings for released and conditional data. Default behavior is to set vintage
    # of data to today's date.
    vint = Dates.format(now(), DSGE_DATE_FORMAT)
    m <= Setting(:data_vintage, vint, true, "vint", "Vintage of data")
    m <= Setting(:cond_vintage, vint, "Vintage of conditional data")
    m <= Setting(:cond_id, "0000", "Identifier of conditional dataset")
    m <= Setting(:cond_full_names, [:obs_gdp, :obs_corepce, :obs_spread, :obs_nominalrate], "Observables used in conditional forecasts")
    m <= Setting(:cond_semi_names, [:obs_spread, :obs_nominalrate], "Observables used in semiconditional forecasts")
    m <= Setting(:use_population_forecast, false, "Whether to use population forecasts as data")
    m <= Setting(:population_mnemonic, :CNP16OV, "Mnemonic of FRED data series for computing per-capita values")

    # Dates
    m <= Setting(:date_presample_start, quartertodate("1959-Q3"), "Start date of pre-sample")
    m <= Setting(:date_prezlb_start, quartertodate("1960-Q1"), "Start date of main sample")
    m <= Setting(:date_zlb_start, quartertodate("2008-Q4"), "Start date of zero lower bound regime")
    m <= Setting(:date_forecast_start, Dates.lastdayofquarter(Dates.today()), "Start date of forecast period")
    m <= Setting(:date_conditional_end, Dates.lastdayofquarter(Dates.today()), "End date of conditional data period")
    m <= Setting(:date_forecast_end, Dates.lastdayofquarter(Dates.today()+Dates.Month(60*3)), "End date of forecast period")

    # Anticipated shocks
    m <= Setting(:n_anticipated_shocks, 0, "Number of anticipated policy shocks")
    m <= Setting(:n_anticipated_shocks_padding, 20, "Padding for anticipated policy shocks")

    # General computation
    m <= Setting(:use_parallel_workers, true, "Use available parallel workers in computations")

    # Estimation
    m <= Setting(:reoptimize, true, "Optimize the posterior mode. If false, reads in mode from a file.")
    m <= Setting(:calculate_hessian, true, "Calculate the hessian at the mode")
    m <= Setting(:n_hessian_test_params, typemax(Int), "Max number of free params for which to calculate Hessian")

    # Metropolis-Hastings
    m <= Setting(:n_mh_simulations, 5000, "Number of draws per block in Metropolis-Hastings")
    m <= Setting(:n_mh_blocks, 22, "Number of blocks for Metropolis-Hastings")
    m <= Setting(:n_mh_burn, 2, "Number of blocks to use as burn-in in Metropolis-Hastings")
    m <= Setting(:mh_thin, 5, "Metropolis-Hastings thinning step")

    # Forecast
    m <= Setting(:forecast_observables, :all, "Observables to forecast")
    m <= Setting(:forecast_kill_shocks, false, "Kill (set to 0) all shocks in forecast")
    m <= Setting(:forecast_tdist_shocks, false, "Draw Students-t distributed shocks in forecast")
    m <= Setting(:forecast_tdist_draw_df, false, "Draw Students-t degrees of freedom parameter")
    m <= Setting(:forecast_tdist_df_val, 15, "Students-t degrees of freedom fixed value")
    m <= Setting(:forecast_smoother, :durbin_koopman, "Choice of smoother to use during forecasting. Can be :kalman, :durbin_koopman, or eventually :carter_kohn")
    m <= Setting(:forecast_jstep, 5, "Forecast thinning step (in addition to MH thinning step")
    m <= Setting(:forecast_enforce_zlb, true, "Enforce zero lower bound in forecast periods")
    m <= Setting(:shockdec_startindex, 190, "Index of start of shock decomposition output period")
    m <= Setting(:shockdec_endindex, 50000, "Index of end of shock decomposition output period")
    m <= Setting(:shockdec_whichshocks, :all, "Sets of shocks for which to conduct shock decomposition")
end

"""
```
default_test_settings!(m::AbstractModel)
```

The following Settings are constructed, initialized and added to
`m.test_settings`. Their purposes are identical to those in
`m.settings`, but these values are used to test DSGE.jl.

### I/O Locations and identifiers
- `saveroot::Setting{ASCIIString}`: A temporary directory in /tmp/
- `dataroot::Setting{ASCIIString}`: dsgeroot/test/reference/
- `data_vintage::Setting{ASCIIString}`: \"_REF\"

### Metropolis-Hastings
- `n_mh_simulations::Setting{Int}`: 100
- `n_mh_blocks::Setting{Int}`: 1
- `n_mh_burn::Setting{Int}`: 0
- `mh_thin::Setting{Int}`: 1
"""
function default_test_settings!(m::AbstractModel)

    test = m.test_settings

    # I/O
    dataroot = normpath(joinpath(dirname(@__FILE__), "..","test","reference"))
    saveroot = mktempdir()

    #General
    test[:saveroot] = Setting(:saveroot, saveroot,
        "Where to write files when in test mode")
    test[:dataroot] = Setting(:dataroot, dataroot,
        "Location of input files when in test mode" )
    test[:data_vintage] = Setting(:data_vintage, "REF", true, "vint",
        "Reference data identifier")
    test[:cond_vintage] = Setting(:cond_vintage, "REF",
        "Vintage of conditional data")
    test[:use_parallel_workers] = Setting(:use_parallel_workers, false, false, "parw",
        "Use available parallel workers in computations")
    test[:n_hessian_test_params] = Setting(:n_hessian_test_params, 3, false, "mhfp",
        "Max number of free params for which to calculate Hessian")

    # Metropolis-Hastings
    test[:n_mh_simulations] = Setting(:n_mh_simulations, 100, false, "nsim",
        "Number of parameter draws per block for testing Metropolis-Hastings")
    test[:n_mh_blocks] = Setting(:n_mh_blocks, 1, false, "nblc",
        "Number of blocks to draw parameters for testing Metropolis-Hastings")
    test[:n_mh_burn] = Setting(:n_mh_burn, 0, false, "nbrn",
        "Number of burn-in blocks for testing Metropolis-Hastings")
    test[:mh_thin] = Setting(:mh_thin, 1, false, "thin",
        "Thinning step for testing Metropolis-Hastings")

    # Forecast
    test[:date_forecast_start] = Setting(:date_forecast_start, quartertodate("2015-Q4"),
        "Start date of forecast period")
    test[:date_forecast_end] = Setting(:date_forecast_end, quartertodate("2016-Q1"),
        "End date of forecast period")
    test[:forecast_jstep] = Setting(:forecast_jstep, 1,
        "Forecast thinning step (in addition to MH thinning step")
    test[:shockdec_startindex] = Setting(:shockdec_startindex, 2,
        "Index of start of shock decomposition output period")
    test[:shockdec_endindex] = Setting(:shockdec_endindex, 4,
        "Index of end of shock decomposition output period")
    test[:shockdec_whichshocks] = Setting(:shockdec_whichshocks, :all, #TODO
        "Sets of shocks for which to conduct shock decomposition")

    return test
end
