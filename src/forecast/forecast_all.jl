"""
```
forecast_all(m::AbstractModel, data::Matrix{Float64}; cond_types::Vector{Symbol}
input_types::Vector{Symbol} output_types::Vector{Symbol}
```

Inputs
------

- `m`: model object
- `data`: matrix of data for observables
- `cond_types`: conditional data type, any combination of
    - `:none`: no conditional data
    - `:semi`: use "semiconditional data" - average of quarter-to-date observations for high frequency series
    - `:full`: use "conditional data" - semiconditional plus nowcasts for desired
      observables
- `input_types`: which set of parameters to use, any combination of
    - `:mode`: forecast using the modal parameters only
    - `:mean`: forecast using the mean parameters only
    - `:full`: forecast using all parameters (full distribution)
    - `:subset`: forecast using a well-defined user-specified subset of draws
- `output_types`: forecast routine outputs to compute, any combination of
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

Outputs
-------

- todo
"""

function forecast_all(m::AbstractModel, data::Matrix{Float64};
                      cond_types::Vector{Symbol}   = Vector{Symbol}(),
                      input_types::Vector{Symbol}  = Vector{Symbol}(),
                      output_types::Vector{Symbol} = Vector{Symbol}())

    # Prepare forecast settings
    use_expected_rate_date     = (n_anticipated_shocks(m) > 0)    # use data on expected future interest rates
    shockdec_whichshocks       = 1:n_shocks_exogenous(m)          # defaults to 1:n_shocks_exogenous(m)
    bounded_interest_rate_flag = true                             # defaults to true, i.e. enforce bounded interest rate
    no_forecast_shocks_flag    = false                            # defaults to false, i.e. do *not* turn off all shocks
    simulation_smoother_flag   = true                             # defaults to true, i.e. do use simulation smoother
    tdist_shocks_flag          = false                            # defaults to false, i.e. do not use t-distributed shocks
    tdist_drawdf_flag          = false                            # defaults to false, i.e. do not draw t-distribution d.o.f.
    tdist_df_value             = 15                               # defaults to 15
    forecast_horizons          = 60                               # defaults to 60
    shockdec_startindex        = 190                              # defaults to index of 2007-Q1, confirm this index is correct
    shockdec_endindex          = size(data,1)+forecast_horizons-1 # defaults to first_forecast_quarter+forecast_horizons-1
    forecast_observables       = m.observables                    # defaults to all observables

    # # Matlab settings to translate
    # ShockremoveList      # indices of shocks to use in shock decomposition and counterfactual
    # bdd_int_rate         # enforce zero lower bound on nominal interest rate in forecast
    # sflag                # set all shocks to zero in forecast
    # dsflag               # draw states from simulation smoother
    # tflag                # draw Shocks from multivariatei t-distribution
    # dfflag               # draw degrees of freedom of t-distributino
    # df_bar               # if !dfflag, set degrees of freedom to value of df_bar (defaults to 15)
    # qahead               # number of periods to forecast past stime
    # Startdate            # index for shock decomosition and counterfactual start date
    # Enddate_forecastfile # total number of quarters in the observable time series and the forecast period

    for input_type in input_types
        for cond_type in cond_types
            for output_type in output_types
                forecast_one(m, data; cond_type=cond_type, input_type=input_type, output_type=output_type)
            end
        end
    end

end

function forecast_one(m::AbstractModel, data::Matrix{Float64};
                      cond_type::Symbol   = :none
                      input_type::Symbol  = :mode
                      output_type::Symbol = :simple)

    # Set up infiles
    input_file_name = get_input_file(m, input_type)

    # Read infiles
    if input_data in [:mean, :mode]
        params = h5open(input_file_name, "r") do f
            read(f, "params")
        end
        TTT = RRR = CCC = zend = []
    elseif input_data == :full
        h5open(input_file_name, "r") do f
            params = read(f, "params")
            TTT = read(f, "TTT")
            RRR = read(f, "RRR")
            CCC = read(f, "CCC")
            zend = read(f, "zend")
        end
    end

    # If we just have one draw of parameters in mode or mean case, then we don't have the
    # pre-computed system matrices. We now recompute them here by running the Kalman filter.
    if input_data in [:mean, :mode]
        update!(m, params)
        sys = compute_system(m; use_expected_rate_date=false)
        kal = filter(m, data, sys; Ny0 = n_presample_periods(m))
        zend = kal[:zend]
    end

    # Prepare conditional data matrix. All missing columns will be set to NaN.
    if cond_type in [:semi, :full]
        cond_data = load_cond_data(m, cond_type)
        data = [data; cond_data]
    end

    # Example: call forecast, unconditional data, states+observables
    forecast(m, data, sys, zend;
             observables=forecast_observables,
             forecast_horizons=forecast_horizons,
             bounded_interest_rate_flag=bounded_interest_rate_flag,
             no_forecast_shocks_flag=no_forecast_shocks_flag,
             tdist_shocks_flag=tdist_shocks_flag,
             tdist_drawdf_flag=tdist_drawdf_flag,
             tdist_df_value=tdist_df_value)

    # Set up outfiles
    output_file_names = get_output_files(m, input_type, output_type, cond_type,
                                              forecast_settings)

    # Write outfiles
    for output_file in output_file_names
        write(output_file)
    end

end

function get_input_file(m, input_type)
    if input_type == :mode
        return rawpath(m,"estimate","paramsmode.h5")
    elseif input_type == :mean
        return workpath(m,"estimate","paramsmean.h5")
    elseif input_type == :full
        return rawpath(m,"estimate","mhsave.h5")
    elseif input_type == :subset
        #TODO
        return ""
    else
        throw(ArgumentError("Invalid input_type: $(input_type)"))
    end
end

function get_output_files(m, input_type, output_type, cond, forecast_settings)
    additional_file_strings = Dict{Symbol,AbstractString}()

    # Add additional file strings here
    additional_file_strings[:para] = abbrev_symbol(input_type)
    additional_file_strings[:cond] = abbrev_symbol(cond)


    filename = rawpath(m, "forecast", output_type * ".h5", additional_file_strings)

    if output_type == :states
        prefix = "histstates"
    elseif output_type == :shocks
        prefix = "histshocks"
    elseif output_type == :shocks_nonstandardized
        prefix = "histshocksns"
    elseif output_type == :forecast
        prefix = "forecaststates"
        prefix = "forecastobs"
        prefix = "forecastshocks"
    elseif output_type == :shockdec
        prefix = "shockdecstates"
        prefix = "shockdecobs"
    elseif output_type == :dettrend
        prefix = "dettrendstates"
        prefix = "dettrendobs"
    elseif output_type == :counter
        prefix = "counterstates"
        prefix = "counterobs"
    elseif output_type in [:simple, :simple_cond]
        prefix = "histstates"
        prefix = "forecaststates"
        prefix = "forecastobs"
        prefix = "forecastshocks"
    elseif output_type == :all

    end
end
