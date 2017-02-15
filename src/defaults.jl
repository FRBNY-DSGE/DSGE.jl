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
    m <= Setting(:use_population_forecast, false, "Whether to use population forecasts as data")
    m <= Setting(:population_mnemonic, :CNP16OV, "Mnemonic of FRED data series for computing per-capita values")

    # Dates
    m <= Setting(:date_presample_start, quartertodate("1959-Q3"), "Start date of pre-sample")
    m <= Setting(:date_mainsample_start, quartertodate("1960-Q1"), "Start date of main sample")
    m <= Setting(:date_zlbregime_start, quartertodate("2008-Q4"), "Start date of zero lower bound regime")
    m <= Setting(:date_mainsample_end, Dates.lastdayofquarter(Dates.today()-Dates.Month(3)), "End date of main sample")
    m <= Setting(:date_forecast_start, Dates.lastdayofquarter(Dates.today()), "Start date of forecast period")
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

    # Sequential Monte Carlo
    m <= Setting(:n_particles, 10000, "Number of particles for use in SMC")
    m <= Setting(:n_Φ, 200, "Number of stages in the tempering schedule")
    m <= Setting(:λ, 3.0, "The 'bending coefficient' λ in Φ(n) = (n/N(Φ))^λ")
    m <= Setting(:n_smc_blocks, 1, "The number of parameter blocks in SMC")
    m <= Setting(:c, .5, "The scaling factor for the covariance, constructed from acpt and tgt during mutation")
    m <= Setting(:init_accept, .25, "The initial average acceptance rate for new particles during mutation")
    m <= Setting(:target_accept, .25, "The initial target acceptance rate for new particles during mutation")

end
