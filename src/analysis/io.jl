################################################################
## I/O to create MeansBands objects
################################################################

"""
```
get_meansbands_input_files(m, input_type, cond_type, output_vars;
                                    [forecast_string = ""], [fileformat = :jld])
```

```
get_meansbands_input_files(directory, filestring_base, input_type, cond_type, output_vars;
                                    [forecast_string = ""], [fileformat = :jld])
```

Returns a dictionary of raw forecast output files to read in to
compute means and bands.

### Input Arguments

#### Method 1:
- `m::AbstractModel`: Model object

#### Method 2:
- `directory::AbstractString`: directory location of input files to read
- `filestring_base::Vector{AbstractString}`: a vector of strings to be
  added as a suffix. These usually come from model settings for which
  print=true. It should *not* include entries for `cond_type` and
  `input_type` (these will be added automatically).

#### Both methods:
- `input_type::Symbol`: See `?forecast_one`
- `cond_type::Symbol`: See `?forecast_one`
- `output_vars::Symbol`: See `?forecast_one`
- `forecast_string::AbstractString`: See `?forecast_one`
- `fileformat`: file extension of saved files
"""
function get_meansbands_input_files(m::AbstractModel, input_type::Symbol,
                                    cond_type::Symbol, output_vars::Vector{Symbol};
                                    forecast_string::AbstractString = "", fileformat = :jld)

    directory = rawpath(m, "forecast")
    base = filestring_base(m)
    get_meansbands_input_files(directory, base, input_type, cond_type, output_vars;
                               forecast_string = forecast_string,
                               fileformat = fileformat)
end
function get_meansbands_input_files(directory::AbstractString, filestring_base::Vector{ASCIIString},
                                    input_type::Symbol, cond_type::Symbol, output_vars::Vector{Symbol};
                                    forecast_string::AbstractString = "", fileformat::Symbol = :jld)

    input_files = Dict{Symbol, ASCIIString}()
    for var in output_vars
        input_files[var] = get_forecast_filename(directory, filestring_base,
                                                 input_type, cond_type, var,
                                                 forecast_string = forecast_string,
                                                 fileformat = fileformat)
        if contains(string(var), "4q")
            input_files[var] = replace(input_files[var], "forecast4q", "forecast")
        end
    end

    return input_files
end

################################################################
## I/O for actual MeansBands objects
################################################################

"""
```
get_meansbands_input_files(m, input_type, cond_type, output_vars;
                                    [forecast_string = ""], [fileformat = :jld])
```

```
get_meansbands_output_files(directory, filestring_base, input_type, cond_type, output_vars;
                                    [forecast_string = ""], [fileformat = :jld])
```

Returns a dictionary of raw forecast output files in which to save
computed means and bands.

### Input Arguments

#### Method 1:
- `m::AbstractModel`: Model object

#### Method 2:
- `directory::AbstractString`: directory location of input files to read
- `filestring_base::Vector{AbstractString}`: a vector of strings to be
  added as a suffix. These usually come from model settings for which
  print=true. It should *not* include entries for `cond_type` and
  `input_type` (these will be added automatically).

#### Both methods:
- `input_type::Symbol`: See `?forecast_one`
- `cond_type::Symbol`: See `?forecast_one`
- `output_vars::Symbol`: See `?forecast_one`
- `forecast_string::AbstractString`: See `?forecast_one`
- `fileformat`: file extension of saved files
"""
function get_meansbands_output_files(m::AbstractModel, input_type::Symbol,
                                     cond_type::Symbol, output_vars::Vector{Symbol};
                                     forecast_string::AbstractString = "",
                                     fileformat::Symbol = :jld)

    directory = workpath(m, "forecast")
    base = filestring_base(m)
    get_meansbands_output_files(directory, base, input_type, cond_type, output_vars;
                                forecast_string = forecast_string, fileformat = fileformat)
end
function get_meansbands_output_files(directory::AbstractString,
                                     filestring_base::Vector{ASCIIString},
                                     input_type::Symbol, cond_type::Symbol, output_vars::Vector{Symbol};
                                     forecast_string::AbstractString = "", fileformat = :jld)

    mb_output_vars = [symbol("mb$x") for x in output_vars]
    output_files = Dict{Symbol, ASCIIString}()
    for var in output_vars
        output_files[var] = get_forecast_filename(directory, filestring_base,
                                                  input_type, cond_type, symbol("mb$var");
                                                  forecast_string = forecast_string,
                                                  fileformat = fileformat)
    end
    return output_files
end


"""
```
function read_mb{S<:AbstractString}(fn::S)
```
Read in a `MeansBands` object saved in `fn`.
"""
function read_mb{S<:AbstractString}(fn::S)
    @assert isfile(fn) "File $fn could not be found"
    mb = jldopen(fn, "r") do f
        read(f, "mb")
    end
    mb
end

################################################################
## I/O for human-readable csvs
################################################################
"""
```
write_meansbands_tables(dirname, mb, [vars = []], [forecast_string = ""],
           [mb_trend = MeansBands()], [mb_dettrend = MeansBands()], [shocks = []])
```

```
write_meansbands_tables(m, mb, [vars = vars], [forecast_string = forecast_string],
              [mb_trend = MeansBands()], [mb_dettrend = MeansBands()], [shocks = []])
```

Write means and bands (ordered) to a csv file in `tablespath(m, "forecast")`

### Inputs
- `m`: Model object
- `mb`: a meansbands object

### Keyword arguments
- `vars`: a `Symbol` or `Vector{Symbol}` indicating which variable in
  `mb` you want to print means and bands for. If left out or empty,
  means and bands for all variables in `mb` will be printed.
- `forecast_string::AbstractString`: See `?forecast_one`
- `mb_trend::MeansBands`: a `MeansBands` object for a trend
  product. This is required when `mb` contains means and bands for
  a shock decomposition.
- `mb_dettrend::MeansBands`: a `MeansBands` object for a deterministic trend
  product. This is required when `mb` contains means and bands for
  a shock decomposition.
- `shocks::Vector{Symbol}`: If `mb` is a shock decomposition, this is
  an optional list of shocks to print to the table. If omitted, all
  shocks will be printed.
"""
function write_meansbands_tables(m::AbstractModel, mb::MeansBands, var::Symbol;
                                 forecast_string::AbstractString = "",
                                 mb_trend::MeansBands = MeansBands(),
                                 mb_dettrend::MeansBands = MeansBands(),
                                 shocks::Vector{Symbol} = Vector{Symbol}())

    # Use model object to compute output directory and settings suffix
    dir = tablespath(m, "forecast")
    base = DSGE.filestring_base(m)

    # Call lower-level function
    write_meansbands_tables(dir, mb, var, shocks = shocks, filestring_base = base,
                            forecast_string = forecast_string,
                            mb_trend = mb_trend, mb_dettrend = mb_dettrend)
end
function write_meansbands_tables(m::AbstractModel, mb::MeansBands;
                                 vars::Vector{Symbol} = Vector{Symbol}(),
                                 mb_trend::MeansBands = MeansBands(),
                                 mb_dettrend::MeansBands = MeansBands(),
                                 shocks::Vector{Symbol} = Vector{Symbol}())

    # Use all vars by default
    if isempty(vars)
        vars = setdiff(names(mb.means), [:date])
    end

    # Write table for each var
    for var in vars
        write_meansbands_tables(m, mb, var,
                                mb_trend = mb_trend, mb_dettrend = mb_dettrend, shocks = shocks)
    end
end
function write_meansbands_tables{S<:AbstractString}(dirname::S, mb::MeansBands;
                                                    vars::Vector{Symbol} = Vector{Symbol}(),
                                                    filestring_base::Vector{S} = Vector{S}(),
                                                    forecast_string::AbstractString = "",
                                                    mb_trend::MeansBands = MeansBands(),
                                                    mb_dettrend::MeansBands = MeansBands(),
                                                    shocks::Vector{Symbol} = Vector{Symbol}())
    # Use all vars by default
    if isempty(vars)
        vars = setdiff(names(mb.means), [:date])
    end

    # Write table for each var
    for var in vars
        write_meansbands_tables(dirname, mb, var, filestring_base = filestring_base,
                                forecast_string = forecast_string,
                                mb_trend = mb_trend, mb_dettrend = mb_dettrend,
                                shocks = shocks)
    end
end
function write_meansbands_tables{S<:AbstractString}(dirname::S, mb::MeansBands, var::Symbol;
                                                    filestring_base::Vector{S} = Vector{S}(),
                                                    forecast_string::AbstractString = "",
                                                    mb_trend::MeansBands = MeansBands(),
                                                    mb_dettrend::MeansBands = MeansBands(),
                                                    shocks::Vector{Symbol} = Vector{Symbol}())
    # What product are we making here?
    prod = get_product(mb)

    # Compute output filename
    filename = "$prod" * "$var.csv"
    filestring_addl = DSGE.get_forecast_filestring_addl(get_para(mb), get_cond_type(mb), forecast_string = forecast_string)
    fullfilename = DSGE.savepath(dirname, filename, filestring_base, filestring_addl)

    # Extract dataframe
    df = if prod in [:hist, :forecast, :forecast4q, :bddforecast, :bddforecast4q]
        prepare_meansbands_table_timeseries(mb, var)
    elseif prod in [:shockdec]
        @assert !isempty(mb_trend)    "Please pass in mb_trend"
        @assert !isempty(mb_dettrend) "Please pass in mb_dettrend"

        prepare_means_table_shockdec(mb, mb_trend, mb_dettrend, var, shocks = shocks)
    end

    # Write to file
    writetable(fullfilename, df)
    println(" * Wrote means and bands for $var to $fullfilename")

    return df
end

"""
```
write_meansbands_tables_all(m, input_type, cond_type, output_vars;
                              [forecast_string = ""], [vars = []], [shocks = shocks])
```

Writes all `output_vars` corresponding to model `m` to tables in `tablespath(m, \"forecast\")`.

### Input arguments
- `m::AbstractModel`: Model object
- `input_type::Symbol`: See `?forecast_one`
- `cond_type::Symbol`: See `?forecast_one`
- `output_vars::Symbol`: See `?forecast_one`

### Keyword arguments
- `forecast_string::AbstractString`: See `?forecast_one`
- `vars`::Vector{Symbol}: Vector of economic variables for which to
  print `output_vars` to tables. If omitted, all shocks will be
  printed.
- `shocks::Vector{Symbol}`: Vector of shocks to print if `output_vars`
  contains a shock decomposition. If omitted, all shocks will be
  printed.
"""
function write_meansbands_tables_all(m::AbstractModel, input_type::Symbol, cond_type::Symbol,
                                     output_vars::Vector{Symbol};
                                     forecast_string = "",
                                     vars::Vector{Symbol} = Vector{Symbol}(),
                                     shocks::Vector{Symbol} = Vector{Symbol}())

    # Create a dictionary of all means and bands
    mbs = Dict{Symbol, MeansBands}()

    # # Figure out stuff for filename
    # base = DSGE.filestring_base(m)
    filestring_addl = DSGE.get_forecast_filestring_addl(input_type, cond_type,
                                                        forecast_string = forecast_string)
    # dirname = tablespath(m, "forecast")
    # filename = DSGE.savepath(dirname, filename, base, filestring_addl)

    # Read in appropriate means and bands
    my_output_vars = add_requisite_output_vars(output_vars)
    for class in [:pseudo, :obs]
        if !isempty(intersect(my_output_vars, [symbol("hist$class"), symbol("histforecast$class")]))
            mbs[symbol("hist$class")] =
                read_mb(workpath(m, "forecast", "mbhist$class.jld", filestring_addl))
        end
        if !isempty(intersect(my_output_vars, [symbol("forecast$class"), symbol("histforecast$class")]))
            mbs[symbol("forecast$class")] =
                read_mb(workpath(m, "forecast", "mbforecast$class.jld", filestring_addl))
        end
        if !isempty(intersect(my_output_vars, [symbol("bddforecast$class"), symbol("bddhistforecast$class")]))
            mbs[symbol("bddforecast$class")] =
                read_mb(workpath(m, "forecast", "mbbddforecast$class.jld", filestring_addl))
        end
        if !isempty(intersect(my_output_vars, [symbol("forecast4q$class"), symbol("histforecast4q$class")]))
            mbs[symbol("forecast4q$class")] =
                read_mb(workpath(m, "forecast", "mbforecast4q$class.jld", filestring_addl))
        end
        if !isempty(intersect(my_output_vars, [symbol("bddforecast4q$class"), symbol("bddhistforecast4q$class")]))
            mbs[symbol("bddforecast4q$class")] =
                read_mb(workpath(m, "forecast", "mbbddforecast4q$class.jld", filestring_addl))
        end
        if symbol("shockdec$class") in my_output_vars
            mbs[symbol("shockdec$class")] =
                read_mb(workpath(m, "forecast", "mbshockdec$class.jld", filestring_addl))
        end
        if symbol("dettrend$class") in my_output_vars
            mbs[symbol("dettrend$class")] =
                read_mb(workpath(m, "forecast", "mbdettrend$class.jld", filestring_addl))
        end
        if symbol("trend$class") in my_output_vars
            mbs[symbol("trend$class")] =
                read_mb(workpath(m, "forecast", "mbtrend$class.jld", filestring_addl))
        end
    end

    # write each output to CSV
    for output in output_vars

        class = get_class(output)

        df = if output in [:histpseudo, :histobs, :forecastpseudo, :forecastobs]

            write_meansbands_tables(m, mbs[output], vars = vars)

        elseif output in [:histforecastpseudo, :histforecastobs]

            mb_histforecast = cat(mb[symbol("hist$(class)")], mb[symbol("forecast$class")])
            write_meansbands_tables(m, mb_histforecast, vars = vars)

        elseif output in [:shockdecpseudo, :shockdecobs]

            mb_shockdec = mbs[symbol("shockdec$(class)")]
            mb_trend    = mbs[symbol("trend$(class)")]
            mb_dettrend = mbs[symbol("dettrend$(class)")]

            write_meansbands_tables(m, mb_shockdec, mb_trend = mb_trend, mb_dettrend = mb_dettrend,
                                    vars = vars, shocks = shocks)
        end
    end
end
