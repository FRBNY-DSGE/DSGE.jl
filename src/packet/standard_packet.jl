"""
```
write_standard_model_packet(m, input_type, cond_type,
    output_vars = [:forecastobs, :forecastpseudo, :shockdecobs, :shockdecpseudo];
    sections = [:estimation, :forecast]
    forecast_string = "", outdir = joinpath(saveroot(m), \"Packets\", spec(m), subspec(m)))
```

Write standard estimation and forecast result packet to `outdir`.
"""
function write_standard_model_packet(m::AbstractDSGEModel, input_type::Symbol, cond_type::Symbol,
                                     output_vars::Vector{Symbol} = [:forecastobs, :forecastpseudo,
                                                                    :shockdecobs, :shockdecpseudo];
                                     sections::Vector{Symbol} = [:estimation, :forecast],
                                     forecast_string::String = "",
                                     outdir::String = joinpath(saveroot(m), "Packets", spec(m), subspec(m)),
                                     purpose::String = "")
    @assert issubset(sections, [:estimation, :forecast, :irf]) "Section specified in `section` kwarg is not supported. Must be a subset of [:estimation, :forecast, :irf]."

    # Title and authors
    title = "Estimation and Forecasting Results. \\\\ " * description(m)
    authors = "User"

    # Compute file name
    cdvt_str = cond_type == :none ? "" : "_cdvt=" * cond_vintage(m)
    fn = joinpath(outdir, "results_cond=" * string(cond_type) * cdvt_str * "_vint=" * data_vintage(m) * ".tex")
    isdir(dirname(fn)) || mkpath(dirname(fn))

    # Write packet
    open(fn, "w") do fid
        write_preamble(fid, title, authors)
        write_spec_section(fid, m, purpose = purpose)
        if :estimation in sections
            write_estimation_section(fid, m)
        end
        if :forecast in sections
            write_forecast_section(fid, m, input_type, cond_type, setdiff(output_vars, [:irfstates, :irfobs, :irfpseudo]),
                                   forecast_string = forecast_string)
        end
        if :irf in sections
            write_irf_section(fid, m, input_type, cond_type,
                              output_vars, forecast_string = forecast_string)
        end
        write_postamble(fid)
    end
    println("Wrote " * fn)
end

"""
```
plot_standard_model_packet(m, input_type, cond_type,
    output_vars = [:forecastobs, :forecastpseudo, :shockdecobs, :shockdecpseudo];
    sections = [:estimation, :forecast]
    forecast_string = "", hist_start-date = 0001-01-01)
```

Plot parameter prior/posterior histograms to `figurespath(m, \"estimate\")` and
forecasts/shock decompositions to `figurespath(m, \"forecast\')`.
"""
function plot_standard_model_packet(m::AbstractDSGEModel, input_type::Symbol, cond_type::Symbol,
                                    output_vars::Vector{Symbol} = [:forecastobs, :forecastpseudo,
                                                                   :shockdecobs, :shockdecpseudo];
                                    sections::Vector{Symbol} = [:estimation, :forecast],
                                    forecast_string::String = "",
                                    hist_start_date::Date = Date("0001-01-01", "yyyy-mm-dd"))
    @assert issubset(sections, [:estimation, :forecast, :irf]) "Section specified in `section` kwarg is not supported. Must be a subset of [:estimation, :forecast, :irf]."
    if :estimation in sections
        plot_prior_posterior(m)
    end
    if :forecast in sections
        for output_var in setdiff(output_vars, [:irfstates, :irfobs, :irfpseudo])
            make_forecast_plots(m, input_type, cond_type, output_var,
                                forecast_string = forecast_string,
                                hist_start_date = hist_start_date)
        end
    end
    if :irf in sections
        plot_irf_section(m, input_type, cond_type, output_vars;
                         forecast_string = forecast_string)
    end
end

"""
```
write_spec_section(fid, m; purpose = "")
```

Write user, model, and data specification section.
"""
function write_spec_section(fid::IOStream, m::AbstractDSGEModel; purpose::String = "")

    @printf fid "\n\n"
    @printf fid "\\section{Specification}\n"
    @printf fid "\n"

    @printf fid "\\subsection{User}\n"
    @printf fid "\n"
    @printf fid "\\begin{itemize}\n"
    @printf fid "  \\item Packet creator: %s\n" splitdir(homedir())[end]
    @printf fid "  \\item Creation time: %s\n" Dates.format(now(), "U d, YYYY at H:MM")
    if !isempty(purpose)
        @printf fid "\\item Purpose: %s\n" purpose
    end
    @printf fid "\\end{itemize}\n"
    @printf fid "\n"

    @printf fid "\\subsection{Model}\n"
    @printf fid "\n"
    @printf fid "\\begin{itemize}\n"
    @printf fid "  \\item Model: %s\n" m.spec
    @printf fid "  \\item Subspec: %s\n" m.subspec
    @printf fid "  \\item Number of states: %d\n" n_states_augmented(m)
    @printf fid "  \\item Number of shocks: %d\n" n_shocks_exogenous(m)
    @printf fid "  \\item Number of observables: %d\n" n_observables(m)
    @printf fid "  \\item Number of parameters: %d\n" n_parameters(m)
    @printf fid "  \\item Number of anticipated shocks: %d\n" n_mon_anticipated_shocks(m)
    @printf fid "\\end{itemize}\n"
    @printf fid "\n"

    @printf fid "\\subsection{Data}\n"
    @printf fid "\n"
    @printf fid "\\begin{itemize}\n"
    @printf fid "  \\item Input data directory: \\path{%s}\n" replace(get_setting(m, :dataroot),"\\"=>"/")
    @printf fid "  \\item Output data directory: \\path{%s}\n" replace(get_setting(m, :saveroot),"\\"=>"/")
    @printf fid "  \\item Data vintage: %s\n" get_setting(m, :data_vintage)
    @printf fid "  \\item Dataset ID: %d\n" get_setting(m, :data_id)
    @printf fid "  \\item Conditional data vintage: %s\n" get_setting(m, :cond_vintage)
    @printf fid "  \\item Conditional dataset ID: %d\n" get_setting(m, :cond_id)
    @printf fid "  \\item Use population forecast: %s\n" get_setting(m, :use_population_forecast)
    @printf fid "  \\item Last quarter of data: %d-Q%d\n" Dates.year(date_mainsample_end(m)) Dates.quarterofyear(date_mainsample_end(m))
    @printf fid "  \\item Last quarter of conditional data: %d-Q%d\n" Dates.year(date_conditional_end(m)) Dates.quarterofyear(date_conditional_end(m))
    @printf fid "\\end{itemize}\n"
end

"""
```
write_estimation_section(fid, m; plotroot = "")
```

Write parameter moment tables and prior/posterior plots.
"""
function write_estimation_section(fid::IOStream, m::AbstractDSGEModel;
                                  plotroot::String = "")
    @printf fid "\n\n"
    @printf fid "\\clearpage\n"
    @printf fid "\\section{Estimation}\n"
    @printf fid "\n"

    @printf fid "\\subsection{Histograms}\n"
    write_estimation_plots(fid, m, plotroot = plotroot)
    @printf fid "\n"

    @printf fid "\\clearpage\n"
    @printf fid "\\subsection{Moments}\n"
    @printf fid "\n"

    @printf fid "\\input{%s}\n" replace(tablespath(m, "estimate", "moments.tex"),"\\"=>"/")
end

"""
```
write_estimation_plots(fid, m; plotroot = "")
```

 Write LaTeX code displaying plots of `m` `product`s to the `IOStream` `fid`. If
`plotroot` is not specified, plots from `figurespath(m, \"estimate\")` are used.
"""
function write_estimation_plots(fid::IOStream, m::AbstractDSGEModel;
                                plotroot::String = "")
    if isempty(plotroot)
        plotroot = replace(figurespath(m, "estimate"), "\\" => "/")
    end
    base = DSGE.filestring_base(m)
    filestr = DSGE.filestring(base, String[])

    @printf fid "\n"
    @printf fid "\\begin{longtable}{cc}\n"
    left = true
    for param in m.parameters
        if param.fixed
            continue
        else
            @printf fid "\\includegraphics[width=0.45\\textwidth]{%s/prior_posterior_%s%s.pdf} %s\n" plotroot DSGE.detexify(param.key) filestr (left ? "&" : "\\\\")
            left = !left
        end
    end

    @printf fid "\\end{longtable}\n"
end

"""
```
write_forecast_section(fid, m, input_type, cond_type,
    output_vars = [:forecastobs, :forecastpseudo, :shockdecobs, :shockdecpseudo];
    forecast_string = "". plotroot = "")
```

Write forecast specification and plots for the given `output_vars`. If `plotroot`
is not specified, plots from `figurespath(m, \"forecast\")` are used.
"""
function write_forecast_section(fid::IOStream, m::AbstractDSGEModel, input_type::Symbol, cond_type::Symbol,
                                output_vars::Vector{Symbol} = [:forecastobs, :forecastpseudo,
                                                               :shockdecobs, :shockdecpseudo];
                                forecast_string::String = "",
                                plotroot::String = "")
    @printf fid "\n\n"
    @printf fid "\\clearpage\n"
    @printf fid "\\section{Forecast}\n"
    @printf fid "\n"

    @printf fid "\\subsection{Specification}\n"
    @printf fid "\n"
    @printf fid "\\begin{itemize}\n"
    @printf fid "  \\item Forecast type: %s\n" input_type
    @printf fid "  \\item Forecast conditional type: %s\n" cond_type
    if !isempty(forecast_string)
        fc_string = if occursin("_", forecast_string)
            replace(forecast_string, "_" => "\\_")
        else
            forecast_string
        end
        @printf fid "  \\item Forecast identifying string: %s\n" fc_string
    end
    @printf fid "  \\item Estimation used: \\path{%s}\n" replace(DSGE.get_forecast_input_file(m, input_type),"\\"=>"/")
    @printf fid "\\end{itemize}\n"

    products = map(get_product, output_vars)

    if :forecast in products
        @printf fid "\\clearpage\n"
        @printf fid "\\subsection{Forecasts}\n"

        for output_var in Base.filter(x -> get_product(x) == :forecast, output_vars)
            write_forecast_plots(fid, m, input_type, cond_type, output_var,
                                 forecast_string = forecast_string,
                                 plotroot = plotroot)
            @printf fid "\n"
        end
    end

    n_output_vars = 2 * (:forecastobs in output_vars) + 4 * (:forecastpseudo in output_vars)
    if n_output_vars > 0
        @printf fid "\\clearpage\n"
        @printf fid "\\subsection{Mean Forecast Estimates}\n"
        write_forecast_table(fid, m, cond_type, output_vars; forecast_string = forecast_string)
    end


    if :shockdec in products
        @printf fid "\\clearpage\n"
        @printf fid "\\subsection{Shock Decompositions}\n"

        for output_var in Base.filter(x -> get_product(x) == :shockdec, output_vars)
            write_forecast_plots(fid, m, input_type, cond_type, output_var,
                                 forecast_string = forecast_string,
                                 plotroot = plotroot)
            @printf fid "\n"
        end
    end
end

"""
```
make_forecast_plots(m, input_type, cond_type, output_var;
    forecast_string = "", plotroot = "", hist_start_date = 0001-01-01,
    trend_nostates = DateFrame())
```

Generate all `output_var` plots for the forecast of `m` specified by the
`input_type`, `cond_type`, and optional `forecast_string`. If `plotroot` is not
specified, plots are saved to `figurespath(m, \"forecast\")`.

trend_nostates is a DataFrame for plotting the trend as a shock in the
shock decomposition. Get the DataFrame via prepare_means_table_trend_nostates.
"""
function make_forecast_plots(m::AbstractDSGEModel, input_type::Symbol, cond_type::Symbol, output_var::Symbol;
                             forecast_string::String = "",
                             plotroot::String = "",
                             hist_start_date::Date = Date("0001-01-01", "yyyy-mm-dd"),
                             trend_nostates::DataFrame = DataFrame(),
                             legend = :bottomright, df_enddate = Date(2100,12,31))

    # Output directory
    if isempty(plotroot)
        plotroot = replace(figurespath(m, "forecast"), "\\" => "/")
    end

    # Set class-specific variables
    class = get_class(output_var)
    product = get_product(output_var)

    if class == :obs
        varnames = collect(keys(m.observables))
        hist_label = "Data"
    elseif class == :pseudo
        varnames = collect(keys(m.pseudo_observables))
        hist_label = "Mean Estimate"
    elseif class == :states
        varnames = collect(keys(m.endogenous_states))
        hist_label = "Mean Estimate"
    else
        error("Unsupported class: " * string(class))
    end

    # Added min() so forecasts that are starting from a year before 2007 will plot properly
    if hist_start_date != Date("0001-01-01", "yyyy-mm-dd")
        start_date = hist_start_date
    else
        start_date = min(DSGE.quartertodate("2007-Q1"), date_mainsample_end(m) - Dates.Year(5))
    end

    if haskey(m.settings, :date_forecast_end)
        end_date = get_setting(m, :date_forecast_end)

        # check forecast horizon was long enough
        end_date -= Dates.Month(3 * get_setting(m, :forecast_horizons))
        if end_date > get_setting(m, :date_forecast_start)
            error("forecast horizons is not long enough for the requested forecast end date")
        else
            end_date = get_setting(m, :date_forecast_end)
        end
    elseif haskey(m.settings, :forecast_horizons)
        end_date = get_setting(m, :date_forecast_start) # this will become the end date
        end_date += Dates.Month((get_setting(m, :forecast_horizons)-1) * 3)
    else
        end_date = DSGE.quartertodate("2020-Q1")
    end

    forecast_label = "Mean Forecast"
    tick_size = 2

    # Make product-specific plot
    product = get_product(output_var)

    if product == :forecast
        plot_history_and_forecast(m, varnames, class, input_type, cond_type,
                                  forecast_string = forecast_string,
                                  bdd_and_unbdd = true,
                                  plotroot = plotroot,
                                  start_date = start_date,
                                  names = Dict(:hist => hist_label, :forecast => forecast_label),
                                  end_date = end_date,
                                  tick_size = tick_size,
                                  legend = legend)
    elseif product == :bddforecast
        plot_history_and_forecast(m, varnames, class, input_type, cond_type,
                                  forecast_string = forecast_string,
                                  bdd_and_unbdd = false, bdd_and_bdd = true,
                                  plotroot = plotroot,
                                  start_date = start_date,
                                  names = Dict(:hist => hist_label, :forecast => forecast_label),
                                  end_date = end_date,
                                  tick_size = tick_size,
                                  legend = legend)
    elseif product == :shockdec
        groups = DSGE.shock_groupings(m)
        plot_shock_decomposition(m, varnames, class, input_type, cond_type,
                                 forecast_string = forecast_string,
                                 plotroot = plotroot,
                                 start_date = start_date,
                                 end_date = end_date,
                                 hist_label = "",
                                 forecast_label = "",
                                 tick_size = tick_size,
                                 legend = legend,
                                 trend_nostates = trend_nostates,
                                 df_enddate = Date(2100,12,31))
    else
        error("Unsupported product: " * string(product))
    end
end

"""
```
write_forecast_plots(fid, m, input_type, cond_type, output_var;
    forecast_string = "", plotroot = "")
```

Write LaTeX code displaying plots of `m` `product`s to the `IOStream` `fid`. If
`plotroot` is not specified, plots from `figurespath(m, \"forecast\")` are used.
"""
function write_forecast_plots(fid::IOStream, m::AbstractDSGEModel,
                              input_type::Symbol, cond_type::Symbol, output_var::Symbol;
                              forecast_string::String = "",
                              plotroot::String = "")
    if isempty(plotroot)
        plotroot = replace(figurespath(m, "forecast"), "\\" => "/")
    end
    base = DSGE.filestring_base(m)
    addl = DSGE.get_forecast_filestring_addl(input_type, cond_type, forecast_string = forecast_string)
    filestr = DSGE.filestring(base, addl)

    class = get_class(output_var)
    if class == :obs
        section_title = "Observables"
        dict = m.observables
    elseif class == :pseudo
        section_title = "Pseudo-Observables"
        dict = m.pseudo_observables
    else
        error("Unsupported class: " * string(class))
    end

    product = get_product(output_var)

    @printf fid "\n"
    @printf fid "\\subsubsection{%s}\n" section_title
    @printf fid "\\begin{longtable}{c}\n"
    for var in keys(dict)
        @printf fid "\\includegraphics[height=0.4\\textheight]{%s/%s_%s%s.pdf} \\\\\n" plotroot product DSGE.detexify(var) filestr
    end
    @printf fid "\\end{longtable}\n"
end

"""
```
write_forecast_table(fid, m, cond_type, output_vars;
    forecast_string = "")
```

Write the standard packet's forecast table listing forecasts 3 years out.

### Inputs

- `fid::IOStream`: file to write to
- `m::AbstractDSGEModel`: model object for current Snapshot forecast
- `cond_type::Symbol`: type of conditional forecast
- `output_vars`: which type of forecast variables to print out

### Keyword Arguments

- `forecast_str::String`: string to use
"""
function write_forecast_table(fid::IOStream, m::AbstractDSGEModel, cond_type::Symbol,
                              output_vars::Vector{Symbol};
                              forecast_string::String = "")

    # Infer dates
    first_fcast_year = Dates.year(date_forecast_start(m))
    fcast_years = [first_fcast_year, first_fcast_year + 1, first_fcast_year + 2, first_fcast_year + 3]
    fcast_dates = map(year -> Date(year, 12, 31), fcast_years)

    curr_month = month_label(m)

    if curr_month != "Jan" && curr_month != "Feb"
        curr_year = "20"*data_vintage(m)[1:2]
    else
        # Because the first forecast of the year still forecasts the previous year.
        curr_year = "20"*string(parse(data_vintage(m)[1:2])-1)
    end

    # Open output file, write header
    year2 = string(Meta.parse(curr_year)+1)
    year3 = string(Meta.parse(curr_year)+2)
    year4 = string(Meta.parse(curr_year)+3)

    write(fid, "\\begin{table}[h!]\n")
    write(fid, "\\centering\n")
    write(fid, "\\begin{tabular}{l | c | c | c | c }\n")
    write(fid, "& " * "$(curr_year)" * " & " * "$(year2)" * " & " * "$(year3)" * " & " * "$(year4)" * " \\\\\n")
    write(fid, "\\hline\n\\hline\n")

    n_output_vars = 2 * (:forecastobs in output_vars) + 4 * (:forecastpseudo in output_vars)
    if :forecastobs in output_vars
        write(fid, print_variable_means(m, cond_type,
                                        :forecast4qobs, :obs_gdp, ["GDP growth ", "(Q4/Q4)"],
                                        fcast_dates, false, forecast_string = forecast_string))
        last_row = n_output_vars == 2 ? true : false
        write(fid, print_variable_means(m, cond_type,
                                        :forecast4qobs, :obs_corepce,
                                        ["Core PCE inflation ", "(Q4/Q4)"],
                                        fcast_dates, last_row, forecast_string = forecast_string))
    end
    if :forecastpseudo in output_vars
        write(fid, print_variable_means(m, cond_type,
                                        :forecast4qpseudo, :π_t, ["Inflation ",
                                                                  "(Q4)"],
                                        fcast_dates, false, forecast_string = forecast_string))
        write(fid, print_variable_means(m, cond_type,
                                        :forecast4qpseudo, :LongRunInflation,
                                        ["Long Run Inflation ", "(Q4)"],
                                        fcast_dates, false, forecast_string = forecast_string))
        write(fid, print_variable_means(m, cond_type,
                                        :forecast4qpseudo, :NaturalRate,
                                        ["Real natural rate", "of interest (Q4)"],
                                        fcast_dates, false, forecast_string = forecast_string))
        write(fid, print_variable_means(m, cond_type,
                                        :forecast4qpseudo, :OutputGap,
                                        ["Output gap", "(Q4)"],
                                        fcast_dates, true, forecast_string = forecast_string))
    end

    write(fid, "\\caption{Mean Forecast Estimates}\n")
    write(fid, "\\end{table}")
end

function print_variable_means(m::AbstractDSGEModel, cond_type::Symbol, output_var::Symbol,
                              varname::Symbol, vardesc::Vector{String}, fcast_dates::Vector{Date},
                              last_row::Bool; forecast_string::String = "")
    n_fcast_dates = length(fcast_dates)

    mb_curr = read_mb(m, :full, cond_type, output_var; forecast_string = forecast_string)

    # Index out Q4 values for current forecast
    all_dates   = mb_curr.metadata[:date_inds]

    # The following works because all_dates is an ordered dictionary so can just treat keys as a array and take their indices.
    # The commented out part stopped working because Julia7 doesn't like indexing with a Date key
    table_dates = findall(x -> x in fcast_dates, collect(keys(all_dates))) #map(date -> all_dates[date], fcast_dates)
    cur_val     = round.(mb_curr.means[:, [:date, varname]][table_dates, 2], digits = 1)

    row1_string = vardesc[1]*repeat(" ", 34-length(vardesc[1])) * "& "
    for i in 1:n_fcast_dates
        if i != n_fcast_dates
            row1_string = string(row1_string, cur_val[i]," & ")
        else
            row1_string = string(row1_string, cur_val[i]," \\\\\n")
        end
    end
    if last_row
        row2_string = vardesc[2]*repeat(" ", 34 - length(vardesc[2]))*repeat("&     ", n_fcast_dates)*"\\\\\n\\end{tabular}"
    else
        row2_string = vardesc[2]*repeat(" ", 34 - length(vardesc[2]))*repeat("&     ", n_fcast_dates)*"\\\\\n\\hline\n"
    end

    return row1_string*row2_string
end


"""
```
write_irf_section(fid, m, input_type, cond_type, forecast_string = "", plotroot = "")
```

Write impulse responses. If `plotroot`
is not specified, plots from `figurespath(m, \"forecast\")` are used.
"""
function write_irf_section(fid::IOStream, m::AbstractDSGEModel, input_type::Symbol, cond_type::Symbol,
                           output_vars::Vector{Symbol}; forecast_string::String = "",
                           plotroot::String = "")
    @printf fid "\n\n"
    @printf fid "\\clearpage\n"
    @printf fid "\\section{Impulse Responses}\n"
    @printf fid "\n"

    @printf fid "\\subsection{Plots}\n"
    if :irfstates in output_vars
        write_irf_plots(fid, m, input_type, cond_type, :irfstates, plotroot = plotroot,
                        forecast_string = forecast_string)
    end
    if :irfobs in output_vars
        write_irf_plots(fid, m, input_type, cond_type, :irfobs, plotroot = plotroot,
                        forecast_string = forecast_string)
    end
    if :irfpseudo in output_vars
        write_irf_plots(fid, m, input_type, cond_type, :irfpseudo, plotroot = plotroot,
                        forecast_string = forecast_string)
    end
    @printf fid "\n"
end

"""
```
write_irf_plots(fid, m, input_type, cond_type, output_var; forecast_string = "", plotroot = "")
```

Plot impulse responses. If `plotroot`
is not specified, plots from `figurespath(m, \"forecast\")` are used.
"""
function write_irf_plots(fid::IOStream, m::AbstractDSGEModel, input_type::Symbol, cond_type::Symbol,
                         output_var::Symbol; forecast_string::String = "", plotroot::String = "")
    if isempty(plotroot)
        plotroot = replace(figurespath(m, "forecast"), "\\" => "/")
    end
    base = DSGE.filestring_base(m)
    addl = DSGE.get_forecast_filestring_addl(input_type, cond_type,
                                             forecast_string = forecast_string)
    filestr = DSGE.filestring(base, addl)

    # Find shocks to remove from list of shocks to compute IRFs
    ant_sh  = [Symbol("rm_shl" * string(i)) for i = 1:n_mon_anticipated_shocks(m)]
    zero_shocks = Vector{Symbol}(undef,0)
    for k in setdiff(keys(m.exogenous_shocks),
                     [Symbol("rm_shl" * string(i)) for i = 1:n_mon_anticipated_shocks(m)])
        if k == :rm_sh && m[:σ_r_m].value == 0. && m[:σ_r_m].fixed
            zero_shocks = vcat(zero_shocks, :rm_sh)
        elseif k == :zp_sh && m[:σ_z_p].value == 0. && m[:σ_z_p].fixed
            zero_shocks = vcat(zero_shocks, :zp_sh)
        elseif k != :rm_sh && k != :zp_sh
            sig_name = Symbol("σ_" * replace(string(k), "_sh" => ""))
            if m[sig_name].value == 0. && m[sig_name].fixed
                zero_shocks = vcat(zero_shocks, k)
            end
        end
    end
    shocks = setdiff(setdiff(collect(keys(m.exogenous_shocks)), ant_sh), zero_shocks)

    class   = get_class(output_var)
    product = get_product(output_var)
    if class == :obs
        section_title = "Observables"
        ant_obs = [Symbol("obs_nominalrate" * string(i)) for i = 1:n_mon_anticipated_shocks(m)]
        rowvars = setdiff(keys(m.observables), ant_obs)
    elseif class == :pseudo
        section_title = "Pseudo-Observables"
        rowvars = keys(m.pseudo_observables)
    elseif class == :states
        section_title = "Endogenous States"
        rowvars = keys(m.endogenous_states)
    else
        error("Unsupported class: " * string(class))
    end

    @printf fid "\n"
    @printf fid "\\subsubsection{%s}\n" section_title
    for var in shocks
        @printf fid "\\begin{figure}[h!]\n"
        @printf fid "\\centering\n"
        @printf fid "\\begin{longtable}{cccc}\n"
        loc = 1 # location within row
        for k in rowvars
            if loc < 4
                endstr = "&"
                loc += 1
            else
                endstr = "\\\\"
                loc = 1
            end
            @printf fid "\\includegraphics[width=0.24\\textwidth]{%s/%s_%s_%s%s.pdf} %s\n" plotroot product DSGE.detexify(var) DSGE.detexify(k) filestr endstr
        end
        @printf fid "\\end{longtable}\n"
        caption_string = if occursin("_", string(var))
            replace(string(DSGE.detexify(var)), "_" => "\\_")
        else
            string(DSGE.detexify(var))
        end
        @printf fid "\\caption{%s}\n" caption_string
        @printf fid "\\end{figure}\n"
    end
    @printf fid "\\clearpage\n"
end

"""
```
plot_irf_section(fid, m, input_type, cond_type; forecast_string = "", plotroot = "", ncols = 4)
```

Plot impulse responses. If `plotroot`
is not specified, plots from `figurespath(m, \"forecast\")` are used.
We assume the default number of columns for the joint plot is 3.
"""
function plot_irf_section(m::AbstractDSGEModel, input_type::Symbol, cond_type::Symbol,
                          output_vars::Vector{Symbol};
                          forecast_string::String = "", plotroot::String = "", ncols::Int64 = 4,
                          kwargs...)
    # Create file name related strings
    if isempty(plotroot)
        plotroot = replace(figurespath(m, "forecast"), "\\" => "/")
    end
    base = DSGE.filestring_base(m)
    addl = DSGE.get_forecast_filestring_addl(input_type, cond_type,
                                             forecast_string = forecast_string)
    filestr = DSGE.filestring(base, addl)


    # Find shocks to remove from list of shocks to compute IRFs
    ant_obs = [Symbol("obs_nominalrate" * string(i)) for i = 1:n_mon_anticipated_shocks(m)]
    ant_sh  = [Symbol("rm_shl" * string(i)) for i = 1:n_mon_anticipated_shocks(m)]
    zero_shocks = Vector{Symbol}(undef,0)
    for k in setdiff(keys(m.exogenous_shocks),
                     [Symbol("rm_shl" * string(i)) for i = 1:n_mon_anticipated_shocks(m)])
        if k == :rm_sh && m[:σ_r_m].value == 0. && m[:σ_r_m].fixed
            zero_shocks = vcat(zero_shocks, :rm_sh)
        elseif k == :zp_sh && m[:σ_z_p].value == 0. && m[:σ_z_p].fixed
            zero_shocks = vcat(zero_shocks, :zp_sh)
        elseif k != :rm_sh && k != :zp_sh
            sig_name = Symbol("σ_" * replace(string(k), "_sh" => ""))
            if m[sig_name].value == 0. && m[sig_name].fixed
                zero_shocks = vcat(zero_shocks, k)
            end
        end
    end

    shocks = setdiff(setdiff(collect(keys(m.exogenous_shocks)), ant_sh), zero_shocks)

    for output_var in intersect(output_vars, [:irfstates, :irfobs, :irfpseudo])
        class = get_class(output_var)
        if class == :obs
            vars = setdiff(collect(keys(m.observables)), ant_obs)
        elseif class == :pseudo
            vars = collect(keys(m.pseudo_observables))
        elseif class == :states
            vars = collect(keys(m.endogenous_states))
        end
        nvars = length(vars)
        nrows = convert(Int, ceil(nvars/ncols))


        for shock in shocks
            # Create individual plots
            plots = plot_impulse_response(m, shock, vars, class, input_type, cond_type,
                                          forecast_string = forecast_string, kwargs...)

            # Plot all variables together
            varplots = if mod(nvars, ncols) > 0
                empty = plot(grid = false, foreground_color_axis = :white, ticks = [])
                n_empty = ncols - mod(nvars, ncols)
                vcat(collect(values(plots)), repeat([empty], n_empty))
            else
                collect(values(plots))
            end
            p = plot(varplots..., layout = grid(nrows, ncols), size = (650, 775))

            # Save
            fn = joinpath(plotroot, "irf_" * string(DSGE.detexify(shock)) * "_" * string(class))
            fn *= filestr * ".pdf"
            DSGE.save_plot(p, fn)
        end
    end
end

# Extends to using two different models
function plot_irf_section(m1::AbstractDSGEModel, m2::AbstractDSGEModel,
                          input_type1::Symbol, input_type2::Symbol,
                          cond_type1::Symbol, cond_type2::Symbol,
                          output_vars::Vector{Symbol};
                          forecast_string1::String = "",
                          forecast_string2::String = "",
                          which_model::Int = 1,
                          plotroot::String = "", ncols::Int64 = 4,
                          bands_pcts::Vector{String} = Vector{String}(undef,0),
                          bands_alpha1::Float64 = 0.5,
                          bands_alpha2::Float64 = 0.5,
                          addl_text::String = "")
    # Create file name related strings
    m = which_model == 1 ? m1 : m2
    if isempty(plotroot)
        plotroot = replace(figurespath(m, "forecast"), "\\" => "/")
    end
    base = DSGE.filestring_base(m)
    addl = DSGE.get_forecast_filestring_addl((which_model == 1) ? input_type1 :
                                             input_type2, (which_model == 2) ?
                                             cond_type1 : cond_type2,
                                             forecast_string = (which_model == 1) ?
                                             forecast_string1 : forecast_string2)
    filestr = DSGE.filestring(base, addl)


    # Find shocks to remove from list of shocks to compute IRFs
    ant_obs = [Symbol("obs_nominalrate" * string(i)) for i = 1:n_mon_anticipated_shocks(m)]
    ant_sh  = [Symbol("rm_shl" * string(i)) for i = 1:n_mon_anticipated_shocks(m)]
    zero_shocks = Vector{Symbol}(undef,0)
    for k in setdiff(keys(m.exogenous_shocks),
                     [Symbol("rm_shl" * string(i)) for i = 1:n_mon_anticipated_shocks(m)])
        if haskey(m.keys, :σ_r_m)
            if k == :rm_sh && m[:σ_r_m].value == 0. && m[:σ_r_m].fixed
                zero_shocks = vcat(zero_shocks, :rm_sh)
            end
        elseif haskey(m.keys, :σ_rm)
            if k == :rm_sh && m[:σ_rm].value == 0. && m[:σ_rm].fixed
                zero_shocks = vcat(zero_shocks, :rm_sh)
            end
        elseif haskey(m.keys, :σ_z_p)
            if k == :zp_sh && m[:σ_z_p].value == 0. && m[:σ_z_p].fixed
                zero_shocks = vcat(zero_shocks, :zp_sh)
            end
        elseif k != :rm_sh && k != :zp_sh
            sig_name = Symbol("σ_" * replace(string(k), "_sh" => ""))
            if m[sig_name].value == 0. && m[sig_name].fixed
                zero_shocks = vcat(zero_shocks, k)
            end
        end
    end

    shocks = setdiff(setdiff(collect(keys(m.exogenous_shocks)), ant_sh), zero_shocks)

    for output_var in intersect(output_vars, [:irfstates, :irfobs, :irfpseudo])
        class = get_class(output_var)
        if class == :obs
            vars = setdiff(collect(keys(m.observables)), ant_obs)
        elseif class == :pseudo
            vars = collect(keys(m.pseudo_observables))
        elseif class == :states
            vars = collect(keys(m.endogenous_states))
        end
        nvars = length(vars)
        nrows = convert(Int, ceil(nvars/ncols))


        for shock in shocks
            # Create individual plots
            plots = plot_impulse_response(m1, m2, shock, vars, class,
                                          input_type1, input_type2,
                                          cond_type1, cond_type2,
                                          forecast_string1 = forecast_string1,
                                          forecast_string2 = forecast_string2,
                                          which_model = which_model,
                                          bands_pcts = bands_pcts,
                                          bands_alpha1 = bands_alpha1,
                                          bands_alpha2 = bands_alpha2,
                                          addl_text = addl_text)

            # Plot all variables together
            varplots = if mod(nvars, ncols) > 0
                empty = plot(grid = false, foreground_color_axis = :white, ticks = [])
                n_empty = ncols - mod(nvars, ncols)
                vcat(collect(values(plots)), repeat([empty], n_empty))
            else
                collect(values(plots))
            end
            p = plot(varplots..., layout = grid(nrows, ncols), size = (650, 775))

            # Save
            fn = joinpath(plotroot, "irf_" * string(DSGE.detexify(shock)) * "_" * string(class) * addl_text)
            fn *= filestr * ".pdf"
            DSGE.save_plot(p, fn)
        end
    end
end
