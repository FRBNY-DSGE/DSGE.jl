"""
```
get_forecast_input_file(m, input_type)
```

Compute the appropriate forecast input filenames for model `m` and
forecast input type `input_type`.

The default input files for each `input_type` can be overriden by adding entries
to the `Dict{Symbol, String}` returned from
`forecast_input_file_overrides(m)`. For example:

```
overrides = forecast_input_file_overrides(m)
overrides[:mode] = \"path/to/input/file.h5\"
```
"""
function get_forecast_input_file(m, input_type;
                                 filestring_addl::Vector{String} = Vector{String}(undef, 0))
    overrides = forecast_input_file_overrides(m)
    if haskey(overrides, input_type)
        override_file = overrides[input_type]
        if ispath(override_file)
            return override_file
        else
            error("Input file $override_file does not exist")
        end
    elseif input_type == :subset
        # If input_type == :subset, also check for existence of overrides[:full]
        return get_forecast_input_file(m, :full)
    end

    if input_type == :mode || input_type == :mode_draw_shocks
        if get_setting(m, :sampling_method) == :MH
            return rawpath(m,"estimate","paramsmode.h5", filestring_addl)
        else
            smc_mode_file = rawpath(m,"estimate","paramsmode.jld2", filestring_addl)
            if isfile(smc_mode_file)
                return smc_mode_file
            else
                return rawpath(m, "estimate", "smc_cloud.jld2", filestring_addl)
            end
        end
    elseif input_type == :mean
        return workpath(m,"estimate","paramsmean.h5", filestring_addl)
    elseif input_type == :init || input_type == :init_draw_shocks
        return ""
    elseif input_type == :prior
        return ""
    elseif input_type in [:full, :subset]
        if get_setting(m, :sampling_method) == :MH
            return rawpath(m, "estimate", "mhsave.h5", filestring_addl)
        elseif get_setting(m, :sampling_method) == :SMC
            return rawpath(m, "estimate", "smcsave.h5", filestring_addl)
        else
            error("Invalid sampling method specification. Change in setting :sampling_method")
        end
    else
        error(ArgumentError("Invalid input_type: $(input_type)"))
    end
end

"""
```
get_forecast_filename(m, input_type, cond_type, output_var;
    pathfcn = rawpath, forecast_string = "", fileformat = :jld2)

get_forecast_filename(directory, filestring_base, input_type, cond_type,
    output_var; forecast_string = "", fileformat = :jld2)
```

### Notes

- If `input_type == :subset`, then the `forecast_string` is also appended to the
filenames. If in this case `forecast_string` is empty, `get_forecast_filename` throws
an error.

- In the second method, `directory` should be a string of the form
  `\"\$saveroot/m990/ss2/forecast/raw/\"`. (Note that a `pathfcn` is therefore
  not required.) `filestring_base` should be equivalent to the result of
  `filestring_base(m)`.
"""
function get_forecast_filename(m::AbstractDSGEModel, input_type::Symbol,
                               cond_type::Symbol, output_var::Symbol;
                               pathfcn::Function = rawpath,
                               forecast_string::String = "",
                               fileformat::Symbol = :jld2)

    # First, we need to make sure we get all of the settings that have been printed to this filestring
    directory = pathfcn(m, "forecast")
    base      = filestring_base(m)

    # FOR IRFS AND OTHER OUTPUT VAR, MAY NEED TO UPDATE SO THAT OUTPUT HAS A REGIME TAG

    # Now we can find the filestring we are looking for
    get_forecast_filename(directory, base, input_type, cond_type, output_var;
                          forecast_string = forecast_string, fileformat = fileformat)
end

function get_forecast_filename(directory::String, filestring_base::Vector{String},
                               input_type::Symbol, cond_type::Symbol, output_var::Symbol;
                               forecast_string::String = "", fileformat = :jld2)

    # Base file name
    file_name = String("$output_var.$fileformat")
    # Gather all of the file strings into an array
    filestring_addl = get_forecast_filestring_addl(input_type, cond_type,
                                                   forecast_string = forecast_string)

    return savepath(directory, file_name, filestring_base, filestring_addl)
end

function get_forecast_filestring_addl(input_type, cond_type; forecast_string::String = "")

    filestring_addl = Vector{String}()
    push!(filestring_addl, String("para=" * String(input_type)))
    push!(filestring_addl, String("cond=" * abbrev_symbol(cond_type)))
    if isempty(forecast_string)
        if input_type == :subset
            error("Must supply a nonempty forecast_string if input_type = subset")
        end
    else
        push!(filestring_addl, String("fcid=" * forecast_string))
    end

    filestring_addl
end

"""
```
get_forecast_output_files(m, input_type, cond_type, output_vars;
    forecast_string = "", fileformat = :jld2)

get_forecast_output_files(directory, filestring_base, input_type, cond_type,
    output_vars; forecast_string = "", fileformat = :jld2)
```

Compute the appropriate output filenames for model `m`, forecast input
type `input_type`, and conditional type `cond_type`, for each output variable in
`output_vars`. Returns a dictionary of file names with one entry for each output_var.

### Arguments

- See `forecast_one` for descriptions of other non-keyword arguments.

### Keyword Arguments

- `forecast_string::String`: subset identifier for when `input_type = :subset`
- `fileformat::Symbol`: file extension, without a period. Defaults to
  `:jld2`, though `:h5` is another common option.

### Notes

See `get_forecast_filename` for more information.
"""
function get_forecast_output_files(m::AbstractDSGEModel, input_type::Symbol,
                                   cond_type::Symbol, output_vars::Vector{Symbol};
                                   forecast_string::String = "",
                                   fileformat::Symbol = :jld2)

    directory = rawpath(m, "forecast")
    base      = filestring_base(m)

    get_forecast_output_files(directory, base, input_type, cond_type, output_vars,
                              forecast_string = forecast_string,
                              fileformat = fileformat)
end

function get_forecast_output_files(directory::String, filestring_base::Vector{String},
                                   input_type::Symbol, cond_type::Symbol, output_vars::Vector{Symbol};
                                   forecast_string::String = "", fileformat::Symbol = :jld2)

    output_files    = Dict{Symbol, String}()
    for var in remove_meansbands_only_output_vars(output_vars)
        output_files[var] = get_forecast_filename(directory, filestring_base,
                                                  input_type, cond_type, var,
                                                  forecast_string = forecast_string,
                                                  fileformat = fileformat)
    end

    return output_files
end

"""
```
write_forecast_outputs(m, output_vars, forecast_output_files, forecast_output;
    df = DataFrame(), block_number = Nullable{Int64}(), block_inds = 1:0,
    verbose = :low)
```

Writes the elements of `forecast_output` indexed by `output_vars` to file, given
`forecast_output_files`, which maps `output_vars` to file names.

If `output_vars` contains `:histobs`, data must be passed in as `df`.
"""
function write_forecast_outputs(m::AbstractDSGEModel, input_type::Symbol,
                                output_vars::Vector{Symbol},
                                forecast_output_files::Dict{Symbol, String},
                                forecast_output::Dict{Symbol, Array{Float64}};
                                df::DataFrame = DataFrame(),
                                block_number::Nullable{Int64} = Nullable{Int64}(),
                                block_inds::AbstractRange{Int64} = 1:0,
                                subset_inds::AbstractRange{Int64} = 1:0,
                                verbose::Symbol = :low)

    for var in output_vars
        prod = get_product(var)
        if prod in [:histut, :hist4q, :forecastut, :bddforecastut, :forecast4q, :bddforecast4q]#, :histlvl, :forecastlvl, :bddforecastlvl]
            # These are computed and saved in means and bands, not
            # during the forecast itself
            continue
        end
        filepath = forecast_output_files[var]
        # Write forecast metadata to a jld2 and the raw forecast output
        # to an h5. The data in the HDF5 will be transferred to the jld2
        # and the h5 file will be deleted when combine_raw_forecast_output_and_metadata
        # is executed.
        if isnull(block_number) || get(block_number) == 1
            JLD2.jldopen(filepath, true, true, true, IOStream) do file
                write_forecast_metadata(m, file, var)
            end
            h5open(replace(filepath, "jld2" => "h5"), "w") do file
                if var == :histobs
                    # :histobs just refers to data, so we only write one draw
                    # (as all draws would be the same)
                    @assert !isempty(df) "df cannot be empty if trying to write :histobs"
                    df1 = df[date_mainsample_start(m) .<= df[!, :date] .<= date_mainsample_end(m), :]
                    data = df_to_matrix(m, df1)

                    # Must call missing2nan since you cannot write Missing types to HDF5 files
                    write(file, "arr", missing2nan(data))
                else
                    # Otherwise, pre-allocate HDF5 dataset which will contain
                    # all draws
                    if !isnull(block_number) && get(block_number) == 1
                        # Determine forecast output size
                        dims = get_forecast_output_dims(m, input_type, var; subset_inds = subset_inds)
                        block_size = forecast_block_size(m)
                        chunk_dims = collect(dims)
                        chunk_dims[1] = block_size

                        # Initialize dataset
                        #pfile = file #.plain
                        dset = isdefined(HDF5, :create_dataset) ?
                            HDF5.create_dataset(file, "arr", datatype(Float64), dataspace(dims...); chunk = chunk_dims) :
                            HDF5.d_create(file, "arr", datatype(Float64), dataspace(dims...), "chunk", chunk_dims)
                    end
                end
            end
        end

        if var != :histobs
            h5open(replace(filepath, "jld2" => "h5"), "r+") do file
                if isnull(block_number)
                    write(file, "arr", Array{Float64}(forecast_output[var]))
                else
                    write_forecast_block(file, Array{Float64}(forecast_output[var]), block_inds)
                end
            end
        end
        println(verbose, :high, " * Wrote $(basename(filepath))")
    end
end

"""
```
write_forecast_metadata(m::AbstractDSGEModel, file::JLDFile, var::Symbol)
```

Write metadata about the saved forecast output `var` to `filepath`.

Specifically, we save dictionaries mapping dates, as well as state, observable,
pseudo-observable, and shock names, to their respective indices in the saved
forecast output array. The saved dictionaries include:

- `date_indices::Dict{Date, Int}`: saved for all forecast outputs except IRFs
- `state_names::Dict{Symbol, Int}`: saved for `var in [:histstates, :forecaststates, :shockdecstates]`
- `observable_names::Dict{Symbol, Int}`: saved for `var in [:forecastobs, :shockdecobs]`
- `observable_revtransforms::Dict{Symbol, Symbol}`: saved identifiers for reverse transforms used for observables
- `pseudoobservable_names::Dict{Symbol, Int}`: saved for `var in [:histpseudo, :forecastpseudo, :shockdecpseudo]`
- `pseudoobservable_revtransforms::Dict{Symbol, Symbol}`: saved identifiers for reverse transforms used for pseudoobservables
- `shock_names::Dict{Symbol, Int}`: saved for `var in [:histshocks, :forecastshocks, :shockdecstates, :shockdecobs, :shockdecpseudo]`

Note that we don't save dates or transformations for impulse response
functions.
"""
function write_forecast_metadata(m::AbstractDSGEModel, file::JLD2.JLDFile, var::Symbol)

    prod  = get_product(var)
    class = get_class(var)

    # Write date range
    if prod != :irf
        dates = if prod == :hist
            quarter_range(date_mainsample_start(m), date_mainsample_end(m))
        elseif prod in [:forecast, :bddforecast]
            quarter_range(date_forecast_start(m), date_forecast_end(m))
        elseif prod in [:shockdec, :dettrend, :trend]
            quarter_range(date_shockdec_start(m), date_shockdec_end(m))
        elseif prod == :decomp
            quarter_range(date_mainsample_start(m), date_forecast_end(m))
        end

        date_indices = Dict(d::Date => i::Int for (i, d) in enumerate(dates))
        write(file, "date_indices", date_indices)
    end

    # Write regimes, if applicable
    if haskey(m.settings, :regime_switching)
        if get_setting(m, :regime_switching)
            write(file, "regime_dates", get_setting(m, :regime_dates)) # Can retrieve regime_indices using date_indices
            if prod == :trend
                write(file, "time_varying_trends", haskey(get_settings(m), :time_varying_trends) ?
                      get_setting(m, :time_varying_trends) : true)
            end
        end
    end

    # Write state names
    if class == :states
        state_indices = merge(m.endogenous_states, m.endogenous_states_augmented)
        @assert length(state_indices) == n_states_augmented(m) # assert no duplicate keys
        write(file, "state_indices", state_indices)
        rev_transforms = Dict{Symbol,Symbol}(x => Symbol("identity") for x in keys(state_indices))
        write(file, "state_revtransforms", rev_transforms)
    end

    # Write observable names and transforms
    if class == :obs
        write(file, "observable_indices", m.observables)
        rev_transforms = if prod != :irf
            Dict{Symbol,Symbol}(x => Symbol(m.observable_mappings[x].rev_transform) for x in keys(m.observables))
        else
            Dict{Symbol,Symbol}(x => Symbol("identity") for x in keys(m.observables))
        end
        write(file, "observable_revtransforms", rev_transforms)
    end

    # Write pseudo-observable names and transforms
    if class == :pseudo
        write(file, "pseudoobservable_indices", m.pseudo_observables)
        rev_transforms = if prod != :irf
            Dict{Symbol,Symbol}(x => Symbol(m.pseudo_observable_mappings[x].rev_transform)
                                for x in keys(m.pseudo_observables))
        else
            Dict{Symbol,Symbol}(x => Symbol("identity") for x in keys(m.pseudo_observables))
        end
        write(file, "pseudoobservable_revtransforms", rev_transforms)
    end

    # Write shock names and transforms
    if class in [:shocks, :stdshocks] || prod in [:shockdec, :irf, :decompshockdec]
        write(file, "shock_indices", m.exogenous_shocks)
        if class in [:shocks, :stdshocks]
            rev_transforms = Dict{Symbol,Symbol}(x => Symbol("identity") for x in keys(m.exogenous_shocks))
            write(file, "shock_revtransforms", rev_transforms)
        end
    end
end

"""
```
write_forecast_block(file::JLDFile, arr::Array, block_number::Int,
    block_inds::AbstractRange{Int64})
```

Writes `arr` to the subarray of `file` indicated by `block_inds`.
"""
function write_forecast_block(file, arr::Array,
                              block_inds::AbstractRange{Int64})
    dataset = isdefined(HDF5, :open_dataset) ? HDF5.open_dataset(file, "arr") : HDF5.d_open(file, "arr")
    dims = size(dataset)
    ndims = length(dims)
    dataset[block_inds, fill(Colon(), ndims-1)...] = arr
    if isdefined(HDF5, :set_extent_dims)
        HDF5.set_extent_dims(dataset, dims)
    else
        HDF5.set_dims!(dataset, dims)
    end
end

"""
```
combine_raw_forecast_output_and_metadata(m, forecast_output_files; verbose = :low)
```
Writes the raw forecast output data (`arr`) saved in the temporary h5 file to the
jld2 file containing the rest of the forecast metadata. The intermediary h5 step exists because
jld2 does not support chunked memory assignment in the same way that jld and h5 permitted previously.
"""
function combine_raw_forecast_output_and_metadata(m::AbstractDSGEModel,
                                                  forecast_output_files::Dict{Symbol, String};
                                                  verbose::Symbol = :low)
    for jld_file_path in values(forecast_output_files)
        h5_file_path = replace(jld_file_path, "jld2" => "h5")
        arr = h5read(h5_file_path, "arr")
        metadata = load(jld_file_path)
        JLD2.jldopen(jld_file_path, true, true, true, IOStream) do file
            write(file, "arr", arr)
            for (k, v) in metadata
                write(file, k, v)
            end
        end
        # Remove the h5 containing the raw forecast output when the data has been transferred.
        rm(h5_file_path)
        if VERBOSITY[verbose] >= VERBOSITY[:low]
            jld_file_name = replace(jld_file_path, rawpath(m, "forecast")*"/" => "")
            h5_file_name = replace(h5_file_path, rawpath(m, "forecast")*"/" => "")
            println("Wrote raw forecast output to $jld_file_name and removed $h5_file_name")
        end
    end
end

"""
```
read_forecast_metadata(file::JLDFile)
```

Read metadata from forecast output files. This includes dictionaries mapping
dates, as well as state, observable, pseudo-observable, and shock names, to
their respective indices in the saved forecast output array. Depending on the
`output_var`, the saved dictionaries might include:

- `date_indices::Dict{Date, Int}`: not saved for IRFs
- `state_indices::Dict{Symbol, Int}`
- `observable_indices::Dict{Symbol, Int}`
- `pseudoobservable_indices::Dict{Symbol, Int}`
- `shock_indices::Dict{Symbol, Int}`
- `state_revtransforms::Dict{Symbol, Symbol}`: states are not transformed, so
  all values are `:identity`
- `observable_revtransforms::Dict{Symbol, Symbol}`
- `pseudoobservable_revtransforms::Dict{Symbol, Symbol}`
- `shock_revtransforms::Dict{Symbol, Symbol}`: shocks are not transformed, so
  all values are `:identity`
"""
function read_forecast_metadata(file::JLD2.JLDFile)
    metadata = Dict{Symbol, Any}()
    for field in setdiff(keys(file), "arr")
        metadata[Symbol(field)] = read(file, field)
    end
    return metadata
end

function read_forecast_output(m::AbstractDSGEModel, input_type::Symbol, cond_type::Symbol,
                              output_var::Symbol, var_name::Symbol,
                              shock_name::Nullable{Symbol} = Nullable{Symbol}();
                              forecast_string::String = "")
    # Get filename
    filename = get_meansbands_input_file(m, input_type, cond_type, output_var,
                                         forecast_string = forecast_string)

    # Determine class and product
    product = get_product(output_var)
    class   = get_class(output_var)

    # Get index corresponding to var_name
    class_long = get_class_longname(class)
    indices = FileIO.load(filename, "$(class_long)_indices")
    var_ind = indices[var_name]

    # Read forecast output
    reg_switch = if haskey(m.settings, :regime_switching)
        get_setting(m, :regime_switching)
    else
        false
    end

    if reg_switch && product == :trend
        fcast_series = read_regime_switching_trend(filename, var_ind)
    else
        fcast_series = if isnull(shock_name)
            read_forecast_series(filename, product, var_ind)
        else
            # Get indices corresponding to shock_name
            shock_name = get(shock_name)
            shock_indices = FileIO.load(filename, "shock_indices")
            shock_ind = shock_indices[shock_name]

            read_forecast_series(filename, var_ind, shock_ind)
        end
    end

    JLD2.jldopen(filename, "r") do file
        # The `fcast_output` for trends only is of size `ndraws`. We
        # need to use `repeat` below because population adjustments will be
        # different in each period. Now we have something of size `ndraws` x `nperiods`
        #
        # If regime switching was used, then `fcast_output` is
        # `ndraws` x `n_regimes`, so we still need to use `repeat`.
        if product == :trend
            if reg_switch
                if !read(file, "time_varying_trends")         # if metadata time_varying_trends is false,
                    regime_dates = read(file, "regime_dates") # then we saved on memory by exploiting fact that
                    date_indices = read(file, "date_indices") # the trends in the state space are not time-varying
                    n_regs       = length(regime_dates)       # but this then requires preprocessing
                    nperiods     = length(date_indices)       # before returning foercast output

                    fcast_series_out = Array{eltype(fcast_series)}(undef, size(fcast_series, 1), nperiods)

                    # Figure out to which regimes the dates in date_indices belong
                    regime_inds = Vector{UnitRange{Int}}(undef, length(regime_dates))
                    trend_dates = sort(collect(keys(date_indices)))
                    for reg in 1:n_regs
                        regime_ind = reg == n_regs ? findall(trend_dates .>= regime_dates[reg]) :
                            findall(regime_dates[reg + 1] .> trend_dates .>= regime_dates[reg])
                        if isnothing(regime_ind) # may be none! so just default to a dummy range
                            regime_inds[reg] = 0:0
                        else
                            regime_inds[reg] = date_indices[trend_dates[regime_ind][1]]:date_indices[trend_dates[regime_ind][end]]
                        end
                    end

                    for (reg, reg_inds) in enumerate(regime_inds)
                        if reg_inds != 0:0
                            fcast_series_out[:, reg_inds] = repeat(fcast_series[:, reg], outer = (1, length(reg_inds)))
                        end
                    end
                    fcast_series = fcast_series_out
                end
            else
                nperiods = length(read(file, "date_indices"))
                fcast_series = repeat(fcast_series, outer = (1, nperiods))
            end
        end

        # Parse transform
        class_long = get_class_longname(class)
        transforms = read(file, string(class_long) * "_revtransforms")
        transform = parse_transform(transforms[var_name])

        # # Handle case when series needs to be accumulated, NEED TO FIX INFLATION B/C FCAST_SERIES MAY BE A MATRIX SOMETIMES, ADD IT INSTEAD AS AN ADDITIONAL DIMENSION -> 3D ARRAY, AND THEN AFTERWARD WE NEED TO ADJUST THE NOMINAL REVERSE TRANSFORMS FOR THIS SYNTAX.
        # if product in [:histlvl, :forecastlvl, :bddforecastlvl]
        #     # Check if inflation is required for accumulation
        #     nominal_accumulation = haskey(get_setting(m, :nominal_accumulated_series), var_name)
        #     if nominal_accumulation
        #         requisite_vars = get_setting(m, var_name) # Must have info on the real series first, then info on the inflation series second
        #         fcast_series = Matrix{eltype(fcast_series)}(undef, length(fcast_series), length(requisite_vars))
        #         for (i, var_info) in enumerate(requisite_vars)
        #             # Get index corresponding to var name
        #             req_var_name, req_var_class = var_info
        #             req_class_long = get_class_longname(req_var_class)
        #             req_indices = FileIO.load(filename, "$(req_class_long)_indices")

        #             # Read forecast series
        #             fcast_series[:, i] = read_forecast_series(filename, product, req_var_ind)
        #         end
        #     end

        #     # Infer appropriate transform for accumulation
        #     transform = get_transformlvl(transform; nominal_transform = nominal_accumulation)
        # end

        fcast_series, transform
    end
end

"""
```
read_forecast_series(filepath, product, var_ind)
read_forecast_series(filepath, var_ind, shock_ind)
```

Read only the forecast output for a particular variable (e.g. for a particular
observable) and possibly a particular shock. Result should be a matrix of size
`ndraws` x `nperiods`.
"""
function read_forecast_series(filepath::String, product::Symbol, var_ind::Int)

    dataset = FileIO.load(filepath, "arr")
    ndims = length(size(dataset))

    # Trends are ndraws x nvars
    if product == :trend
        whole = FileIO.load(filepath, "arr")
        if ndims == 1 # one draw
            arr = whole[var_ind, Colon()]
            arr = reshape(arr, (1, 1))
        elseif ndims == 2 # many draws
            arr = whole[Colon(), var_ind]
            arr = reshape(arr, (length(arr), 1))
        end

    # Other products are ndraws x nvars x nperiods
    elseif product in [:hist, :histut, :hist4q, :forecast, :forecastut, :forecast4q,
                       :bddforecast, :bddforecastut, :bddforecast4q, :dettrend,
                       :decompdata, :decompnews, :decomppara, :decompdettrend, :decomptotal]
                       #:forecastlvl, :histlvl, :bddforecastlvl]
        inds_to_read = if ndims == 2 # one draw
            whole = FileIO.load(filepath, "arr")
            arr = whole[var_ind, Colon()]
            arr = reshape(arr, (1, length(arr)))
        elseif ndims == 3 # many draws
            whole = FileIO.load(filepath, "arr")
            arr = whole[Colon(), var_ind, Colon()]
        end
    else
        error("Invalid product: $product for this method")
    end

    return arr
end

function read_forecast_series(filepath::String, var_ind::Int, shock_ind::Int)

    dataset = FileIO.load(filepath, "arr")
    ndims = length(size(dataset))

    if ndims == 3 # one draw
        whole = FileIO.load(filepath, "arr")
        arr = whole[var_ind, :, shock_ind]
        arr = reshape(arr, 1, length(arr))
    elseif ndims == 4 # many draws
        whole = FileIO.load(filepath, "arr")
        arr = whole[:, var_ind, :, shock_ind]
    end

    return arr
end

"""
```
read_regime_switching_trend(filepath, var_ind)
```

Read only the trend output for a particular variable (e.g. for a particular
observable). Result should be a matrix of size `ndraws` × `n_regimes` or
`ndraws` × `n_data_periods`, depending on whether the state space
system was time-varying or not.
"""
function read_regime_switching_trend(filepath::String, var_ind::Int)
    whole = FileIO.load(filepath, "arr")
    if ndims(whole) == 2 # one draw
        arr = whole[var_ind, Colon()]
        arr = reshape(arr, (1, length(arr)))
    elseif ndims(whole) == 3 # many draws
        arr = whole[Colon(), var_ind, Colon()]
    end
    return arr
end
