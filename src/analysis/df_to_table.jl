"""
`function construct_fcast_and_hist_dfs(m, cond_type, vars, save_to_table, table_caption,
filename, savedir, forecast_string, include_T_in_df_forecast, use_4q, hist_start,
forecast_end)`

Construct and return two DataFrame objects that respectively contain the forecast and corresponding
history series for a set of specified variables. Alternatively, if `save_to_table` is
specified to be true with `table_caption`, `filename`, and `savedir` also specified, then
`construct_fcast_and_hist_dfs` will save the DataFrames as LaTeX tables in the specified directory.

### Arguments
- `m::AbstractModel`: The model object that was used for the forecast.
- `cond_type::Symbol`: The conditioning of the forecast.
- `vars::Vector{Symbol}`: The vector of forecasted variables (can be any combination of observables or pseudo-observables) from which to construct the table.

### Keyword Arguments
- `save_to_table::Bool`: Boolean indicator of whether or not to save the table. If set to be true, then one should also specify the `table_caption`, `filename`, and `savedir`.
- `forecast_string::String`: The forecast string (if relevant) for the given forecast.
- `include_T_in_df_forecast`: Whether or not the last historical period value should be included in the forecast table for reference.
- `use_4q::Bool`: Whether to pull the 4q output_var as opposed to the standard quarterly.
- `hist_start::Date`: The date from which to start the historical table.
- `forecast_end::Date`: The date from which to end the forecast table.
"""
function construct_fcast_and_hist_dfs(m::AbstractModel, cond_type::Symbol,
                                      vars::Vector{Symbol};
                                      save_to_table::Bool = false,
                                      table_caption::String = "",
                                      filename::String = "",
                                      savedir::String = "",
                                      forecast_string::String = "",
                                      include_T_in_df_forecast::Bool = true,
                                      use_4q::Bool = false,
                                      hist_start::Date = DSGE.quartertodate("2008-Q4"),
                                      forecast_end::Date =
                                      DSGE.iterate_quarters(m.settings[:date_forecast_start].value, 12))

    # Assert that if save_to_table is true, there is a caption, filename, and a savedir
    if save_to_table
        @assert !isempty(table_caption) && !isempty(filename) && !isempty(savedir)
    end

    if use_4q
        mb_histobs      = read_mb(m, :full, cond_type, :hist4qobs, forecast_string = forecast_string)
        mb_histpseudo   = read_mb(m, :full, cond_type, :hist4qpseudo, forecast_string = forecast_string)
        mb_forecastobs      = read_mb(m, :full, cond_type, :forecast4qobs, forecast_string = forecast_string)
        mb_forecastpseudo   = read_mb(m, :full, cond_type, :forecast4qpseudo, forecast_string = forecast_string)

        mb_histobs = create_q4q4_mb(mb_histobs)
        mb_histpseudo = create_q4q4_mb(mb_histpseudo)
        mb_forecastobs = create_q4q4_mb(mb_forecastobs)
        mb_forecastpseudo = create_q4q4_mb(mb_forecastpseudo)
    else
        mb_histobs      = read_mb(m, :full, cond_type, :histobs, forecast_string = forecast_string)
        mb_histpseudo   = read_mb(m, :full, cond_type, :histpseudo, forecast_string = forecast_string)
        mb_forecastobs      = read_mb(m, :full, cond_type, :forecastobs, forecast_string = forecast_string)
        mb_forecastpseudo   = read_mb(m, :full, cond_type, :forecastpseudo, forecast_string = forecast_string)
    end

    obs = intersect(vars, m.observables.keys)
    pseudo = intersect(vars, m.pseudo_observables.keys)

    df_histobs    = mb_histobs.means[:, vcat(:date, obs)]
    df_histpseudo = mb_histpseudo.means[:, vcat(:date, pseudo)]
    hist_start_ind = findfirst(x -> x == hist_start, df_histobs[:date])
    df_histobs    = df_histobs[hist_start_ind:end, :]
    df_histpseudo = df_histpseudo[hist_start_ind:end, :]

    df_forecastobs      = mb_forecastobs.means[:, vcat(:date, obs)]
    df_forecastpseudo   = mb_forecastpseudo.means[:, vcat(:date, pseudo)]
    if use_4q
        # If producing 4q figures, then the forecast_end date must be a Q4 date
        forecast_end = DSGE.quartertodate(string(Dates.year(forecast_end))*"-Q4")
    end
    forecast_end_ind = findfirst(x -> x == forecast_end, df_forecastobs[:date])
    df_forecastobs      = df_forecastobs[1:forecast_end_ind, :]
    df_forecastpseudo   = df_forecastpseudo[1:forecast_end_ind, :]

    if include_T_in_df_forecast
        df_forecastobs = append!(df_histobs[end, :], df_forecastobs)
        df_forecastpseudo = append!(df_histpseudo[end, :], df_forecastpseudo)
    end

    df_forecast = join(df_forecastobs, df_forecastpseudo, on = :date)
    df_hist     = join(df_histobs, df_histpseudo, on =:date)

    obs_keys = m.observables.keys
    pseudo_keys = m.pseudo_observables.keys

    # Constructs mapping from standard key names in observable/pseudo observable means-bands objects
    # to actual names of variables and the corresponding inverse mapping (for the purposes
    # of creating a mapping of actual names of variables to units)
    header_mappings = create_table_header_mappings(m, vars)
    units = unit_mappings(m, df_forecast, header_mappings, use_4q = use_4q)

    rename!(df_forecast, header_mappings)
    rename!(df_hist, header_mappings)

    if save_to_table
        fcast_table_caption =  table_caption*" Forecast"
        hist_table_caption  =  table_caption*" History"
        fcast_filename = filename*"_forecast"
        hist_filename  = filename*"_history"

        df_to_table(df_forecast, fcast_table_caption, fcast_filename, savedir, units)
        df_to_table(df_hist, hist_table_caption, hist_filename, savedir, units)
    else
        return df_forecast, df_hist
    end
end

# Implement a df_to_table like function that splits up the df/units into
# sub-DataFrames/Dictionaries of 3/4 variables and creates a single LaTeX
# document that has all of the variables in it

# Enforce that the first column of df is a date column named :date

"""
`function df_to_table(df, caption, filename, savedir, units)`

This is the low level function that is called by `construct_fcast_and_hist_dfs` if the
`save_to_table` kwarg for that function is set to be true.

Alternatively, this function can be called directly, provided with the relevant arguments
for naming, labeling, and saving the table.

### Arguments
- `df::DataFrame`: The DataFrame object that is storing the various series.
- `caption::String`: The title of the LaTeX table.
- `filename::String`: The name of the file.
- `savedir::String`: The filepath ending in the lowest level directory that should contain the table.
- `units::OrderedDict{Symbol, String}`: A dictionary that maps the column names of `df` to
the units of that particular series (this calculation is done automatically in if
`save_to_table` is set to be true.
"""
function df_to_table(df::DataFrame, caption::String, filename::String, savedir::String,
                     units::OrderedDict{Symbol, String})

    # Open the TeX file
    savedir = savedir[end] == "/" ? savedir : savedir*"/"
    if !ispath(savedir)
        println("The savedir path provided does not currently exist. Do you want to create the path '"*savedir*"'? y/n")
        answer = readline(STDIN)
        if answer == "y"
            mkpath(savedir)
        else
            throw("Create the proper savedir and call df_to_table again.")
        end
    else
        mkpath(savedir)
    end
    table_out = savedir*filename*".tex"
    fid = open(table_out, "w")

    DSGE.write_table_preamble(fid)
    function write_single_table(fid::IOStream, df::DataFrame, units::OrderedDict{Symbol, String})

        # Write header
        n_columns = length(df.columns)
        col_str   = repeat("c", n_columns)

        @printf fid "%s%s%s" "\\begin{longtable}{" col_str "}\n"
        @printf fid "\\caption{%s}\n" caption
        @printf fid "\\\\ \\hline\n"

        # Write column names
        date_range = df[:date].data

        column_keys = names(df)

        @printf fid "%s " column_keys[1]
        for i in 2:n_columns
            column_key = string(column_keys[i])
            if ismatch(r"_", column_key) # if the key has an underscore then replace it with the proper LaTeX syntax
                sub_strs = split(column_key, "_")
                column_key = sub_strs[1]*"\\_"*sub_strs[2]
            end
            column_entry = "\\parbox\\{0.3\\linewidth\\}\\{\\centering "*column_key*"\\}"
            if i != n_columns
                @printf fid "& %s " column_entry
            else
                @printf fid "& %s \\\\\n" column_entry
            end
        end

        # Write units
        for i in 2:n_columns
            if i != n_columns
                @printf fid "& %s " units[column_keys[i]]
            else
                @printf fid "& %s\n" units[column_keys[i]]
            end
        end
        @printf fid "\\\\ \\hline\n"
        @printf fid "\\endhead\n"

        for (i, date) in enumerate(date_range)
            for (j, key) in enumerate(column_keys)
                if j != length(column_keys)
                    if key == :date
                        @printf fid "%s " df[key][i]
                    else
                        @printf fid "& %.2f " df[key][i]
                    end
                else
                    @printf fid "& %.2f \\\\\n" df[key][i]
                end
            end
        end

        @printf fid "\\end{longtable}\n"
    end

    k = 1
    for i in 1:3:length(units)
        units_keys   = i+2 < length(units) ? units.keys[i:i+2] : units.keys[i:end]
        df_subset    = DataFrame()
        units_subset = OrderedDict{Symbol, String}()

        df_subset[:date] = df[:date]
        for unit in units_keys
            df_subset[unit] = df[unit]
            units_subset[unit] = units[unit]
        end
        write_single_table(fid, df_subset, units_subset)
        if k % 2 == 0
            @printf fid "\\clearpage\n" # every two tables, break the page
        else
            @printf fid "\\vspace*{.5cm}\n"
        end
        k += 1
    end

    @printf fid "\\end{document}"

    close(fid)
end

# Rename keys in the obs dictionaries
# So that the DataFrame has LaTeX conformant names
function create_table_header_mappings(m::AbstractModel, vars::Vector{Symbol})
    obs_keys = m.observables.keys
    header_mappings = OrderedDict{Symbol, Symbol}()
    for var in vars
        if var in obs_keys
            header_mappings[var] = Symbol(m.observable_mappings[var].name)
        else
            header_mappings[var] = DSGE.detexify(Symbol(m.pseudo_observable_mappings[var].name))
        end
    end
    return header_mappings
end

# Defining the units for the variables included
function unit_mappings(m::AbstractModel, df::DataFrame,
                       header_mappings::OrderedDict{Symbol, Symbol}; use_4q::Bool = false)
    units = OrderedDict{Symbol, String}()
    quarter = use_4q ? "Q4" : "Q"
    obs_keys = m.observables.keys
    pseudo_keys = m.pseudo_observables.keys
    for key in names(df)
        name = key != :date ? header_mappings[key] : continue
        if key in obs_keys
            if ismatch(r"pct_annualized", string(m.observable_mappings[key].rev_transform))
                units[name] = "("*quarter*"/"*quarter*") \\% Annualized"
            elseif ismatch(r"quartertoannual", string(m.observable_mappings[key].rev_transform))
                units[name] = quarter
            end
        elseif key in pseudo_keys
            if ismatch(r"pct_annualized", string(m.pseudo_observable_mappings[key].rev_transform))
                units[name] = "("*quarter*"/"*quarter*") \\% Annualized"
            elseif ismatch(r"quartertoannual", string(m.pseudo_observable_mappings[key].rev_transform))
                units[name] = quarter
            elseif ismatch(r"identity", string(m.pseudo_observable_mappings[key].rev_transform))
                units[name] = quarter
            end
        end
    end
    return units
end
