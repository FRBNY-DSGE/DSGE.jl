################################################################
## I/O to create MeansBands objects
################################################################

"""
```
get_meansbands_input_files(m, input_type, cond_type, output_vars;
    forecast_string = "", fileformat = :jld)
```

```
get_meansbands_input_files(directory, filestring_base, input_type, cond_type, output_vars;
    forecast_string = "", fileformat = :jld)
```

Returns a dictionary of raw forecast output files to read in to compute means
and bands.

### Inputs

**Method 1:**

- `m::AbstractModel`

**Method 2:**

- `directory::String`: directory location of input files to read
- `filestring_base::Vector{String}`: a vector of strings to be added as
  a suffix. These usually come from model settings for which `print = true`. It
  should *not* include entries for `cond_type` and `input_type` (these will be
  added automatically).

**Both methods:**

- `input_type::Symbol`: See `?forecast_one`
- `cond_type::Symbol`: See `?forecast_one`
- `output_vars::Symbol`: See `?forecast_one`
- `forecast_string::String`: See `?forecast_one`
- `fileformat`: file extension of saved files
"""
function get_meansbands_input_files(m::AbstractModel, input_type::Symbol,
                                    cond_type::Symbol, output_vars::Vector{Symbol};
                                    forecast_string::String = "", fileformat = :jld)

    directory = rawpath(m, "forecast")
    base = filestring_base(m)
    get_meansbands_input_files(directory, base, input_type, cond_type, output_vars;
                               forecast_string = forecast_string,
                               fileformat = fileformat)
end

function get_meansbands_input_files(directory::String, filestring_base::Vector{String},
                                    input_type::Symbol, cond_type::Symbol, output_vars::Vector{Symbol};
                                    forecast_string::String = "", fileformat::Symbol = :jld)

    input_files = Dict{Symbol, String}()
    for var in output_vars
        input_files[var] = get_forecast_filename(directory, filestring_base,
                                                 input_type, cond_type, var,
                                                 forecast_string = forecast_string,
                                                 fileformat = fileformat)
        if contains(string(var), "forecast4q")
            input_files[var] = replace(input_files[var], "forecast4q", "forecast")
        elseif contains(string(var), "hist4q")
            input_files[var] = replace(input_files[var], "hist4q", "hist")
        end
    end

    return input_files
end


################################################################
## I/O for actual MeansBands objects
################################################################

"""
```
get_meansbands_output_files(m, input_type, cond_type, output_vars;
    forecast_string = "", fileformat = :jld)
```

```
get_meansbands_output_files(directory, filestring_base, input_type, cond_type, output_vars;
    forecast_string = "", fileformat = :jld)
```

Returns a dictionary of raw forecast output files in which to save
computed means and bands.

### Inputs

**Method 1:**

- `m::AbstractModel`: Model object

**Method 2:**

- `directory::String`: directory location of input files to read
- `filestring_base::Vector{String}`: a vector of strings to be
  added as a suffix. These usually come from model settings for which
  print=true. It should *not* include entries for `cond_type` and
  `input_type` (these will be added automatically).

**Both methods:**

- `input_type::Symbol`: See `?forecast_one`
- `cond_type::Symbol`: See `?forecast_one`
- `output_vars::Symbol`: See `?forecast_one`
- `forecast_string::String`: See `?forecast_one`
- `fileformat`: file extension of saved files
"""
function get_meansbands_output_files(m::AbstractModel, input_type::Symbol,
                                     cond_type::Symbol, output_vars::Vector{Symbol};
                                     forecast_string::String = "",
                                     fileformat::Symbol = :jld)

    directory = workpath(m, "forecast")
    base = filestring_base(m)
    get_meansbands_output_files(directory, base, input_type, cond_type, output_vars;
                                forecast_string = forecast_string, fileformat = fileformat)
end

function get_meansbands_output_files(directory::String,
                                     filestring_base::Vector{String},
                                     input_type::Symbol, cond_type::Symbol, output_vars::Vector{Symbol};
                                     forecast_string::String = "", fileformat = :jld)

    mb_output_vars = [Symbol("mb$x") for x in output_vars]
    output_files = Dict{Symbol, String}()
    for var in output_vars
        output_files[var] = get_forecast_filename(directory, filestring_base,
                                                  input_type, cond_type, Symbol("mb$var");
                                                  forecast_string = forecast_string,
                                                  fileformat = fileformat)
    end
    return output_files
end

"""
```
read_mb(fn::String)

read_mb(m, input_type, cond_type, output_var; forecast_string = "",
    bdd_and_unbdd::Bool = false)
```

Read in a `MeansBands` object saved in `fn`, or use the model object `m` to
determine the file location.

If `bdd_and_unbdd`, then `output_var` must be either `:forecast` or
`:forecast4q`. Then this function calls `read_bdd_and_unbdd` to return a
`MeansBands` with unbounded means and bounded bands.
"""
function read_mb(fn::String)
    @assert isfile(fn) "File $fn could not be found"
    jldopen(fn, "r") do f
        read(f, "mb")
    end
end

function read_mb(m::AbstractModel, input_type::Symbol, cond_type::Symbol,
                 output_var::Symbol; forecast_string::String = "",
                 bdd_and_unbdd::Bool = false)

    if bdd_and_unbdd
        @assert get_product(output_var) in [:forecast, :forecast4q]
        bdd_output_var = Symbol(:bdd, output_var)
        files = get_meansbands_output_files(m, input_type, cond_type,
                                            [output_var, bdd_output_var];
                                            forecast_string = forecast_string)
        read_bdd_and_unbdd_mb(files[bdd_output_var], files[output_var])
    else
        files = get_meansbands_output_files(m, input_type, cond_type, [output_var];
                                            forecast_string = forecast_string)
        read_mb(files[output_var])
    end
end

"""
```
read_bdd_and_unbdd_mb(bdd_fn::String, unbdd_fn::String)
```

Read in the bounded and unbounded forecast `MeansBands` from `bdd_fn` and
`unbdd_fn`. Create and return a `MeansBands` with the unbounded means and
bounded bands.
"""
function read_bdd_and_unbdd_mb(bdd_fn::String, unbdd_fn::String)
    # Check files exist
    @assert isfile(bdd_fn)   "File $bdd_fn could not be found"
    @assert isfile(unbdd_fn) "File $unbdd_fn could not be found"

    # Read MeansBands
    bdd_mb   = read_mb(bdd_fn)
    unbdd_mb = read_mb(unbdd_fn)

    # Check well-formed
    for fld in [:para, :forecast_string, :cond_type, :date_inds, :class, :indices]
        @assert bdd_mb.metadata[fld] == unbdd_mb.metadata[fld] "$fld field does not match: $((bdd_mb.metadata[fld], unbdd_mb.metadata[fld]))"
    end
    @assert (bdd_mb.metadata[:product], unbdd_mb.metadata[:product]) in [(:bddforecast, :forecast), (:bddforecast4q, :forecast4q)] "Invalid product fields: $((bdd_mb.metadata[:product], unbdd_mb.metadata[:product]))"

    # Stick together unbounded means, bounded bands
    return MeansBands(unbdd_mb.metadata, unbdd_mb.means, bdd_mb.bands)
end


################################################################
## I/O for human-readable csvs
################################################################

"""
```
write_meansbands_tables_timeseries(m, input_type, cond_type, output_var;
    forecast_string = "", bdd_and_unbdd = false,
    dirname = tablespath(m, \"forecast\"), kwargs...)

write_meansbands_tables_timeseries(dirname, filestring_base, mb;
    tablevars = get_variables(mb))
```

### Inputs

**Method 1 only:**

- `m::AbstractModel`
- `input_type::Symbol`
- `cond_type::Symbol`
- `output_var::Symbol`: `class(output_var)` must be one of `[:hist, :forecast, :hist4q, :forecast4q, :bddforecast, :bddforecast4q, :trend, :dettrend, :histforecast, :histforecast4q]`

**Method 2 only:**

- `dirname::String`: directory to which tables are saved
- `filestring_base::Vector{String}`: the result of `filestring_base(m)`,
  typically `[\"vint=yymmdd\"]``


### Keyword Arguments

- `tablevars::Vector{Symbol}`: which series to write tables for

**Method 1 only:**

- `forecast_string::String`
- `bdd_and_unbdd::Bool`: whether to use unbounded means and bounded
  bands. Applies only for `class(output_var) in [:forecast, :forecast4q]`
- `dirname::String`: directory to which tables are saved
"""
function write_meansbands_tables_timeseries(m::AbstractModel, input_type::Symbol,
                                            cond_type::Symbol, output_var::Symbol;
                                            forecast_string::String = "",
                                            bdd_and_unbdd::Bool = false,
                                            dirname::String = tablespath(m, "forecast"),
                                            kwargs...)
    # Check that output_var is a time series
    prod = get_product(output_var)
    @assert prod in [:hist, :forecast, :hist4q, :forecast4q, :bddforecast, :bddforecast4q,
                     :trend, :dettrend, :histforecast, :histforecast4q]

    # Read in MeansBands
    class = get_class(output_var)

    if prod in [:histforecast, :histforecast4q]
        fourq = contains(string(prod), "4q")
        mb_hist     = read_mb(m, input_type, cond_type, Symbol(fourq ? :hist4q : :hist, class),
                              forecast_string = forecast_string)
        mb_forecast = read_mb(m, input_type, cond_type, Symbol(fourq ? :forecast4q : :forecast, class),
                              forecast_string = forecast_string, bdd_and_unbdd = bdd_and_unbdd)
        mb = cat(mb_hist, mb_forecast)
    else
        mb = read_mb(m, input_type, cond_type, output_var, forecast_string = forecast_string,
                     bdd_and_unbdd = bdd_and_unbdd)
    end

    # Call second method
    write_meansbands_tables_timeseries(dirname, filestring_base(m), mb;
                                       kwargs...)
end

function write_meansbands_tables_timeseries(dirname::String, filestring_base::Vector{String},
                                            mb::MeansBands;
                                            tablevars::Vector{Symbol} = get_variables(mb))
    for tablevar in tablevars
        df = prepare_meansbands_table_timeseries(mb, tablevar)
        write_meansbands_table(dirname, filestring_base, mb, df, tablevar)
    end
end

"""
```
write_means_tables_shockdec(m, input_type, cond_type, class;
    forecast_string = "", dirname = tablespath(m, \"forecast\"),
    kwargs...)

write_means_tables_shockdec(dirname, filestring_base, mb_shockdec,
    mb_trend, mb_dettrend, mb_hist, mb_forecast; tablevars = get_variables(mb),
    columnvars = get_shocks(mb))
```

### Inputs

**Method 1 only:**

- `m::AbstractModel`
- `input_type::Symbol`
- `cond_type::Symbol`
- `class::Symbol`

**Method 2 only:**

- `dirname::String`: directory to which tables are saved
- `filestring_base::Vector{String}`: the result of `filestring_base(m)`,
  typically `[\"vint=yymmdd\"]``
- `mb_shockdec::MeansBands`
- `mb_trend::MeansBands`
- `mb_dettrend::MeansBands`
- `mb_hist::MeansBands`: optional
- `mb_forecast::MeansBands`: optional

### Keyword Arguments

- `tablevars::Vector{Symbol}`: which series to write tables for
- `columnvars::Vector{Symbol}`: which shocks to include as columns in the tables

**Method 1 only:**

- `forecast_string::String`
- `bdd_and_unbdd::Bool`: whether to use unbounded means and bounded
  bands. Applies only for `class(output_var) in [:forecast, :forecast4q]`
- `dirname::String`: directory to which tables are saved
"""
function write_means_tables_shockdec(m::AbstractModel, input_type::Symbol,
                                     cond_type::Symbol, class::Symbol;
                                     forecast_string::String = "",
                                     dirname::String = tablespath(m, "forecast"),
                                     kwargs...)
    # Read in necessary MeansBands
    products = [:shockdec, :trend, :dettrend, :hist, :forecast]
    output_vars = [Symbol(product, class) for product in products]
    mbs = OrderedDict{Symbol, MeansBands}()
    for output_var in output_vars
        try
            mbs[output_var] = read_mb(m, input_type, cond_type, output_var, forecast_string = forecast_string)
        catch ex
            if output_var in [Symbol(:hist, class), Symbol(:forecast, class)]
                warn("MeansBands for " * string(output_var) * " not found")
                mbs[output_var] = MeansBands()
            else
                rethrow(ex)
            end
        end
    end

    # Call second method
    write_means_tables_shockdec(dirname, filestring_base(m), values(mbs)...;
                                kwargs...)
end

function write_means_tables_shockdec(dirname::String, filestring_base::Vector{String},
                                     mb_shockdec::MeansBands, mb_trend::MeansBands,
                                     mb_dettrend::MeansBands,
                                     mb_hist::MeansBands = MeansBands(),
                                     mb_forecast::MeansBands = MeansBands();
                                     tablevars::Vector{Symbol} = get_variables(mb_shockdec),
                                     columnvars::Vector{Symbol} = get_shocks(mb_shockdec))
    for tablevar in tablevars
        df = prepare_means_table_shockdec(mb_shockdec, mb_trend, mb_dettrend, tablevar,
                                          mb_forecast = mb_forecast, mb_hist = mb_hist,
                                          shocks = columnvars)
        write_meansbands_table(dirname, filestring_base, mb_shockdec, df, tablevar)
    end
end

"""
```
write_meansbands_tables_irf(m, input_type, cond_type, class;
    forecast_string = "", dirname = tablespath(m, \"forecast\"),
    kwargs...)

write_meansbands_tables_irf(dirname, filestring_base, mb;
    tablevars = get_shocks(mb), columnvars = get_variables(mb))
```

### Inputs

**Method 1 only:**

- `m::AbstractModel`
- `input_type::Symbol`
- `cond_type::Symbol`
- `class::Symbol`

**Method 2 only:**

- `dirname::String`: directory to which tables are saved
- `filestring_base::Vector{String}`: the result of `filestring_base(m)`,
  typically `[\"vint=yymmdd\"]``
- `mb::MeansBands`

### Keyword Arguments

- `tablevars::Vector{Symbol}`: which shocks to write tables for
- `columnvars::Vector{Symbol}`: which series' impulse responses to include as
  columns in the tables

**Method 1 only:**

- `forecast_string::String`
- `bdd_and_unbdd::Bool`: whether to use unbounded means and bounded
  bands. Applies only for `class(output_var) in [:forecast, :forecast4q]`
- `dirname::String`: directory to which tables are saved
"""
function write_meansbands_tables_irf(m::AbstractModel, input_type::Symbol,
                                     cond_type::Symbol, class::Symbol;
                                     forecast_string::String = "",
                                     dirname::String = tablespath(m, "forecast"),
                                     kwargs...)
    output_var = Symbol(:irf, class)

    # Read in MeansBands
    mb = read_mb(m, input_type, cond_type, output_var, forecast_string = forecast_string)

    # Call second method
    write_meansbands_tables_irf(dirname, filestring_base(m), mbs...;
                                kwargs...)
end

function write_meansbands_tables_irf(dirname::String, filestring_base::Vector{String},
                                     mb::MeansBands,
                                     tablevars::Vector{Symbol} = get_shocks(mb),
                                     columnvars::Vector{Symbol} = get_variables(mb))
    for tablevar in tablevars
        df = prepare_meansbands_table_irf(mb, tablevar, columnvars)
        write_meansbands_table(dirname, filestring_base, mb, df, tablevar)
    end
end

"""
```
write_meansbands_table(dirname, filestring_base, mb, df, tablevar)
```

### Inputs

- `dirname::String`: directory to which tables are saved. Defaults to
  `tablespath(m, \"forecast\")`
- `filestring_base::Vector{String}`: the result of `filestring_base(m)`,
  typically `[\"vint=yymmdd\"]``
- `mb::MeansBands`: used for computing the output file name
- `df::DataFrame`: the result of calling one of
  `prepare_meansbands_table_timeseries`, `prepare_means_table_shockdec`, or
  `prepare_means_table, irf`
- `tablevar::Symbol`: used for computing the base output file name
"""
function write_meansbands_table(dirname::String, filestring_base::Vector{String},
                                mb::MeansBands, df::DataFrame, tablevar::Symbol)
    # Extract metadata
    prod = get_product(mb)
    para = get_para(mb)
    cond = get_cond_type(mb)
    forecast_string = mb.metadata[:forecast_string]

    # Compute output file name
    filename = detexify(string(prod) * "_" * string(tablevar) * ".csv")
    filestring_addl = get_forecast_filestring_addl(para, cond, forecast_string = forecast_string)
    fullfilename = savepath(dirname, filename, filestring_base, filestring_addl)

    # Write to file
    writetable(fullfilename, df)
    println(" * Wrote means and bands for $tablevar to $fullfilename")
end

"""
```
write_meansbands_tables_all(m, input_type, cond_type, output_vars;
    forecast_string = "", dirname = tablespath(m, \"forecast\"),
    vars = [], shocks = shocks)
```

Write all `output_vars` corresponding to model `m` to tables in `dirname`.

### Inputs

- `m::AbstractModel`
- `input_type::Symbol`: See `?forecast_one`
- `cond_type::Symbol`: See `?forecast_one`
- `output_vars::Symbol`: See `?forecast_one`

### Keyword Arguments

- `forecast_string::String`: See `?forecast_one`
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
                                     bdd_and_unbdd::Bool = false,
                                     dirname::String = tablespath(m, "forecast"),
                                     vars::Vector{Symbol} = Vector{Symbol}(),
                                     shocks::Vector{Symbol} = Vector{Symbol}())
    for output_var in output_vars

        class = get_class(output)
        prod  = get_product(output)

        if prod in [:hist, :forecast, :hist4q, :forecast4q, :bddforecast, :bddforecast4q,
                    :trend, :dettrend, :histforecast, :histforecast4q]
            write_meansbands_table_timeseries(m, input_type, cond_type, output_var,
                                              tablevars = vars, columnvars = shocks,
                                              forecast_string = forecast_string,
                                              bdd_and_unbdd = bdd_and_unbdd,
                                              dirname = dirname)

        elseif prod == :shockdec
            write_means_tables_shockdec(m, input_type, cond_type, class,
                                        tablevars = vars, columnvars = shocks,
                                        bdd_and_unbdd = bdd_and_unbdd,
                                        forecast_string = forecast_string,
                                        dirname = dirname)

        elseif prod == :irf
            write_means_tables(m, input_type, cond_type, class,
                               tablevars = shocks, columnvars = vars,
                               forecast_string = forecast_string,
                               dirname = dirname)
        end
    end
end