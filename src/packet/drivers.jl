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
usual_estimation(m; reoptimize = true, paramsstart_file = "")
```

Load data, estimate, and write parameter moment tables. If `reoptimize`,
optimization can optionally start from a provided `paramsstart_file`. If
`!reoptimize`, the mode and hessian will be read in from `rawpath(m,
\"estimate\")`.
"""
function usual_estimation(m::AbstractModel;
                          reoptimize::Bool = true,
                          paramsstart_file::String = "")
    if reoptimize && !isempty(paramsstart_file)
        params = h5read(paramsstart_file, "params")
        DSGE.update!(m, params)
    else
        specify_mode!(m, get_forecast_input_file(m, :mode))
        specify_hessian(m, rawpath(m, "estimate", "hessian.h5"))
    end

    df = load_data(m)
    estimate(m, df)
    moment_tables(m)
end

"""
```
usual_forecast(m, input_type, cond_type,
    output_vars = [:histobs, :histpseudo, :forecastobs, :forecastpseudo];
    est_override = "", forecast_string = "",
    density_bands = [0.5, 0.6, 0.7, 0.8, 0.9],
    mb_matrix = false)
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
                        shock_var_value::Float64 = 0.0)

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
                 shock_var_value = shock_var_value)

    # Compute means and bands
    compute_meansbands(m, input_type, cond_type, output_vars; forecast_string = forecast_string,
                       density_bands = density_bands, verbose = :high)

    # Convert MeansBands to matrices if desired
    if mb_matrix
        meansbands_to_matrix(m, input_type, cond_type, output_vars; forecast_string = forecast_string,
                             verbose = :high)
    end
end

function altpolicy_settings!(m::AbstractModel, altpolicy::AltPolicy)
    # I/O and data settings
    fn = dirname(@__FILE__)
    dataroot = "$(fn)/save/input_data/"
    saveroot = "$(fn)/save/"
    m <= Setting(:dataroot, dataroot, "Input data directory path")
    m <= Setting(:saveroot, saveroot, "Output data directory path")
    m <= Setting(:use_population_forecast, true)
    m <= Setting(:alternative_policy, altpolicy, true, "apol", "Alternative policy")
end

function make_altpolicy(key::Symbol)
    # Eval-ing a Symbol gets you a function
    eval(key)()
end

function use_estimation_vintage!(m::AbstractModel, input_type::Symbol)
    # Determine estimation vintage
    est_file = get_forecast_input_file(m, input_type)
    est_vint = match(r"\d{6}", basename(est_file)).match

    # Cache current vintage and update to estimation vintage
    fcast_vint = data_vintage(m)
    m <= Setting(:data_vintage, est_vint)
    return fcast_vint
end

function usual_decompose_forecast(m_new::M, m_old::M, df_new::DataFrame, df_old::DataFrame,
                            input_type::Symbol, cond_new::Symbol, cond_old::Symbol,
                            classes::Vector{Symbol}; est_override::String = "",
                            verbose::Symbol = :low, forecast_string_old = "", forecast_string_new = "", kwargs...) where M<:AbstractModel

    if !isempty(est_override)
        overrides = forecast_input_file_overrides(m)
        overrides[input_type] = est_override
        overrides = forecast_input_file_overrides(m_old)
        overrides[input_type] = est_override
    end

    decompose_forecast(m_new, m_old, df_new, df_old, input_type, cond_new, cond_old, classes; verbose = verbose, forecast_string_old = forecast_string_old, forecast_string_new = forecast_string_new)

end
