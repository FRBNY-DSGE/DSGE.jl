"""
```
forecast_all(...)
```

Inputs
------

- `m`: model object
- `data`: matrix of data for observables
- `cond`: conditional data type, any combination of
    - `:none`: no conditional data
    - `:semi`: use "semiconditional data" - average of quarter-to-date observations for high frequency series
    - `:cond`: use "conditional data" - semiconditional plus nowcasts for desired
      observables
- `output_type`: any combination of
    - `:states`: smoothed states (history) for all specified conditional data types
    - `:shocks`: smoothed shocks (history, standardized) for all specified conditional data types
    - `:shocks_nonstandardized`: smoothed shocks (history, non-standardized) for all
        specified conditional data types
    - `:forecast`: forecast of states and observables for all specified conditional data types
    - `:shockdec`: shock decompositions (history) of states and observables for all
        specified conditional data types
    - `:dettrend`: deterministic trend (history) of states and observables for all specified
        conditional data types
    - `:counter`: counterfactuals (history) of states and observables for all specified
        conditional data types
    - `:simple`: smoothed states, forecast of states, forecast of observables for
        *unconditional* data only
    - `:simple_cond`: smoothed states, forecast of states, forecast of observables for all
        specified conditional data types
    - `:all`: smoothed states (history), smoothed shocks (history, standardized), smoothed
      shocks (history, non-standardized), shock decompositions (history), deterministic
      trend (history), counterfactuals (history), forecast, forecast shocks drawn, shock
      decompositions (forecast), deterministic trend (forecast), counterfactuals (forecast)
   Note that some similar outputs may or may not fall under the "forecast_all" framework,
   including
    - `:mats`: recompute system matrices (TTT, RRR, CCC) given parameters only
    - `:zend`: recompute final state vector (s_{T}) given parameters only
    - `:irfs`: impulse response functions
"""

function forecast_all(m::AbstractModel, data::Matrix{Float64};
                      cond::Vector{Symbol}    = Vector{Symbol}(),
                      output_type::Vector{Symbol} = Vector{Symbol}())

    # Prepare forecast settings
    zerobound                  = (n_anticipated_shocks(m) > 0)
    shockdec_whichshocks       = 1:n_shocks_exogenous(m)          # defaults to 1:n_shocks_exogenous(m)
    bounded_interest_rate_flag = 1                                # defaults to 1, i.e. enforce bounded interest rate
    no_forecast_shocks_flag    = 0                                # defaults to 0, i.e. do *not* turn off all shocks
    disturbance_smoother_flag  = 1                                # defaults to 1, i.e. do use disturbance smoother
    tdist_shocks_flag          = 0                                # defaults to 0, i.e. do not use t-distributed shocks
    tdist_drawdf_flag          = 0                                # defaults to 0, i.e. do not draw t-distribution d.o.f.
    tdist_df_value             = 15                               # defaults to 15
    forecast_horizons          = 60                               # defaults to 60
    shockdec_startindex        = 190                              # defaults to index of 2007-Q1, confirm this index is correct
    shockdec_endindex          = size(data,1)+forecast_horizons-1 # defaults to first_forecast_quarter+forecast_horizons-1

    # # Matlab settings to translate
    # ShockremoveList      # indices of shocks to use in shock decomposition and counterfactual
    # bdd_int_rate         # enforce zero lower bound on nominal interest rate in forecast
    # sflag                # set all shocks to zero in forecast
    # dsflag               # draw states from disturbance smoother
    # tflag                # draw Shocks from multivariatei t-distribution
    # dfflag               # draw degrees of freedom of t-distributino
    # df_bar               # if !dfflag, set degrees of freedom to value of df_bar (defaults to 15)
    # qahead               # number of periods to forecast past stime
    # Startdate            # index for shock decomosition and counterfactual start date
    # Enddate_forecastfile # total number of quarters in the observable time series and the forecast period

    # Set up infiles
    input_file_names = set_forecast_infiles(m)

    # Read infiles
    for input_file in input_file_names
        read(input_file)
    end

    # Prepare conditional data matrices
    # May not be necessary as we can iterate over it

    # Call forecast level 2 programs
    for (cond_data_type, cond_data) in cond
        data_for_forecast = append!(copy(data), cond_data)
        forecast(m, data_for_forecast...)
    end

    # Set up outfiles
    output_file_names = set_forecast_outfiles(m, output_type)

    # Write outfiles
    for output_file in output_file_names
        write(output_file)
    end

end
