"""
```
usual_settings!(m, vint; cdvt = vint, dsid = data_id(m), cdid = cond_id(m),
    fcast_date = Dates.lastdayofquarter(Dates.today()),
    altpolicy = AltPolicy(:historical, eqcond, solve))
```

Apply usual defaults for the following settings:

- `data_vintage` and `cond_vintage`: given by input argument `vint`
- `date_forecast_start` and `date_conditional_end`: given by kwarg `fcast_date`
- `use_population_forecast`: `true`
- `alternative_policy`: given by input argument `altpolicy`. If this argument is
  specified, then `altpolicy_settings!` and `altpolicy.setup` are also called.
"""
function usual_settings!(m::AbstractModel, vint::String;
                         cdvt::String = vint,
                         dsid::Int = data_id(m),
                         cdid::Int = cond_id(m),
                         fcast_date::Dates.Date = Dates.lastdayofquarter(Dates.today()),
                         altpolicy::AltPolicy = AltPolicy(:historical, eqcond, solve))
    m <= Setting(:data_vintage, vint)
    m <= Setting(:cond_vintage, cdvt)
    m <= Setting(:data_id, dsid)
    m <= Setting(:cond_id, cdid)
    m <= Setting(:date_forecast_start,  fcast_date)
    m <= Setting(:date_conditional_end, fcast_date)
    m <= Setting(:use_population_forecast, true)

    if altpolicy.key != :historical
        altpolicy_settings!(m, altpolicy)
    end
end

"""
```
usual_forecast(m, input_type, cond_type,
    output_vars = [:histobs, :histpseudo, :forecastobs, :forecastpseudo];
    est_override = "", forecast_string = "",
    density_bands = [0.5, 0.6, 0.7, 0.8, 0.9],
    mb_matrix = false, check_empty_columns = true)
```

Forecast, compute means and bands, and optionally (if `mb_matrix`) convert
`MeansBands` to matrices. If the path `est_override` is provided, it will be
added to `forecast_input_file_overrides(m)`.
"""
function usual_forecast(m::AbstractModel, input_type::Symbol, cond_type::Symbol,
                        output_vars::Vector{Symbol} = [:histobs, :histpseudo, :forecastobs, :forecastpseudo,
                                                       :shockdecobs, :shockdecpseudo];
                        est_override::String = "",
                        forecast_string::String = "",
                        density_bands::Vector{Float64} = [0.5, 0.6, 0.7, 0.8, 0.9],
                        mb_matrix::Bool = false,
                        shock_name::Symbol = :none,
                        shock_var_name::Symbol = :none,
                        shock_var_value::Float64 = 0.0,
                        check_empty_columns = true)

    # Override estimation file if necessary
    if !isempty(est_override)
        overrides = forecast_input_file_overrides(m)
        overrides[input_type] = est_override
    end

    # Forecast
    forecast_one(m, input_type, cond_type, output_vars; forecast_string = forecast_string,
                 verbose = :high,
                 shock_name = shock_name,
                 shock_var_name = shock_var_name,
                 shock_var_value = shock_var_value,
                 check_empty_columns = check_empty_columns)

    # Compute means and bands
    compute_meansbands(m, input_type, cond_type, output_vars; forecast_string = forecast_string,
                       density_bands = density_bands, check_empty_columns = check_empty_columns,
                       verbose = :high)

    # Convert MeansBands to matrices if desired
    if mb_matrix
        meansbands_to_matrix(m, input_type, cond_type, output_vars; forecast_string = forecast_string,
                             verbose = :high)
    end
end
