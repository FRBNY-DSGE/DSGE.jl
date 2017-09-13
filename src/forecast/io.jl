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
function get_forecast_input_file(m, input_type)
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

    if input_type == :mode
        return rawpath(m,"estimate","paramsmode.h5")
    elseif input_type == :mean
        return workpath(m,"estimate","paramsmean.h5")
    elseif input_type == :init
        return ""
    elseif input_type in [:full, :subset]
        return rawpath(m,"estimate","mhsave.h5")
    else
        throw(ArgumentError("Invalid input_type: $(input_type)"))
    end
end

"""
```
get_forecast_filename(m, input_type, cond_type, output_var;
    pathfcn = rawpath, forecast_string = "", fileformat = :jld)

get_forecast_filename(directory, filestring_base, input_type, cond_type,
    output_var; forecast_string = "", fileformat = :jld)
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
function get_forecast_filename(m::AbstractModel, input_type::Symbol,
                               cond_type::Symbol, output_var::Symbol;
                               pathfcn::Function = rawpath,
                               forecast_string::String = "",
                               fileformat::Symbol = :jld)

    # First, we need to make sure we get all of the settings that have been printed to this filestring
    directory = pathfcn(m, "forecast")
    base      = filestring_base(m)

    # Now we can find the filestring we are looking for
    get_forecast_filename(directory, base, input_type, cond_type, output_var;
                          forecast_string = forecast_string, fileformat = fileformat)
end

function get_forecast_filename(directory::String, filestring_base::Vector{String},
                               input_type::Symbol, cond_type::Symbol, output_var::Symbol;
                               forecast_string::String = "", fileformat = :jld)

    # Base file name
    file_name = String("$output_var.$fileformat")

    # Gather all of the file strings into an array
    filestring_addl = get_forecast_filestring_addl(input_type, cond_type,
                                                   forecast_string = forecast_string)

    return savepath(directory, file_name, filestring_base, filestring_addl)
end

function get_forecast_filestring_addl(input_type, cond_type; forecast_string::String = "")

    filestring_addl = Vector{String}()
    push!(filestring_addl, String("para=" * abbrev_symbol(input_type)))
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
    forecast_string = "", fileformat = :jld)

get_forecast_output_files(directory, filestring_base, input_type, cond_type,
    output_vars; forecast_string = "", fileformat = :jld)
```

Compute the appropriate output filenames for model `m`, forecast input
type `input_type`, and conditional type `cond_type`, for each output variable in
`output_vars`. Returns a dictionary of file names with one entry for each output_var.

### Arguments

- See `forecast_one` for descriptions of other non-keyword arguments.

### Keyword Arguments

- `forecast_string::String`: subset identifier for when `input_type = :subset`
- `fileformat::Symbol`: file extension, without a period. Defaults to
  `:jld`, though `:h5` is another common option.

### Notes

See `get_forecast_filename` for more information.
"""
function get_forecast_output_files(m::AbstractModel, input_type::Symbol,
                                   cond_type::Symbol, output_vars::Vector{Symbol};
                                   forecast_string::String = "",
                                   fileformat::Symbol = :jld)

    directory = rawpath(m, "forecast")
    base      = filestring_base(m)

    get_forecast_output_files(directory, base, input_type, cond_type, output_vars,
                              forecast_string = forecast_string,
                              fileformat = fileformat)
end

function get_forecast_output_files(directory::String, filestring_base::Vector{String},
                                   input_type::Symbol, cond_type::Symbol, output_vars::Vector{Symbol};
                                   forecast_string::String = "", fileformat::Symbol = :jld)

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
function write_forecast_outputs(m::AbstractModel, input_type::Symbol,
                                output_vars::Vector{Symbol},
                                forecast_output_files::Dict{Symbol, String},
                                forecast_output::Dict{Symbol, Array{Float64}};
                                df::DataFrame = DataFrame(),
                                block_number::Nullable{Int64} = Nullable{Int64}(),
                                block_inds::Range{Int64} = 1:0,
                                subset_inds::Range{Int64} = 1:0,
                                verbose::Symbol = :low)

    for var in output_vars
        prod = get_product(var)
        if prod in [:hist4q, :forecast4q, :bddforecast4q]
            # These are computed and saved in means and bands, not
            # during the forecast itself
            continue
        end
        filepath = forecast_output_files[var]
        if isnull(block_number) || get(block_number) == 1
            jldopen(filepath, "w") do file
                write_forecast_metadata(m, file, var)

                if var == :histobs
                    # :histobs just refers to data, so we only write one draw
                    # (as all draws would be the same)
                    @assert !isempty(df) "df cannot be empty if trying to write :histobs"
                    data = df_to_matrix(m, df; include_presample = false)
                    write(file, "arr", data)

                else
                    # Otherwise, pre-allocate HDF5 dataset which will contain
                    # all draws
                    if !isnull(block_number) && get(block_number) == 1
                        # Determine forecast output size
                        dims  = get_forecast_output_dims(m, input_type, var; subset_inds = subset_inds)
                        block_size = forecast_block_size(m)
                        chunk_dims = collect(dims)
                        chunk_dims[1] = block_size

                        # Initialize dataset
                        pfile = file.plain
                        HDF5.d_create(pfile, "arr", datatype(Float64), dataspace(dims...), "chunk", chunk_dims)
                    end
                end
            end
        end

        if var != :histobs
            jldopen(filepath, "r+") do file
                if isnull(block_number)
                    write(file, "arr", forecast_output[var])
                else
                    write_forecast_block(file, forecast_output[var], block_inds)
                end
            end
        end

        if VERBOSITY[verbose] >= VERBOSITY[:high]
            println(" * Wrote $(basename(filepath))")
        end
    end
end

"""
```
write_forecast_metadata(m::AbstractModel, file::JldFile, var::Symbol)
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

Note that we don't save dates or transformations for impulse response functions.
"""
function write_forecast_metadata(m::AbstractModel, file::JLD.JldFile, var::Symbol)

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
        end

        date_indices = Dict(d::Date => i::Int for (i, d) in enumerate(dates))
        write(file, "date_indices", date_indices)
    end

    # Write state names
    if class == :states
        state_indices = merge(m.endogenous_states, m.endogenous_states_augmented)
        @assert length(state_indices) == n_states_augmented(m) # assert no duplicate keys
        write(file, "state_indices", state_indices)
        rev_transforms = Dict{Symbol,Symbol}(x => Symbol("DSGE.identity") for x in keys(state_indices))
        write(file, "state_revtransforms", rev_transforms)
    end

    # Write observable names and transforms
    if class == :obs
        write(file, "observable_indices", m.observables)
        rev_transforms = if prod != :irf
            Dict{Symbol,Symbol}(x => Symbol(m.observable_mappings[x].rev_transform) for x in keys(m.observables))
        else
            Dict{Symbol,Symbol}(x => Symbol("DSGE.identity") for x in keys(m.observables))
        end
        write(file, "observable_revtransforms", rev_transforms)
    end

    # Write pseudo-observable names and transforms
    if class == :pseudo
        pseudo, pseudo_mapping = pseudo_measurement(m)
        write(file, "pseudoobservable_indices", pseudo_mapping.inds)
        rev_transforms = if prod != :irf
            Dict{Symbol,Symbol}(x => Symbol(pseudo[x].rev_transform) for x in keys(pseudo))
        else
            Dict{Symbol,Symbol}(x => Symbol("DSGE.identity") for x in keys(pseudo))
        end
        write(file, "pseudoobservable_revtransforms", rev_transforms)
    end

    # Write shock names and transforms
    if class in [:shocks, :stdshocks] || prod in [:shockdec, :irf]
        write(file, "shock_indices", m.exogenous_shocks)
        if class in [:shocks, :stdshocks]
            rev_transforms = Dict{Symbol,Symbol}(x => Symbol("DSGE.identity") for x in keys(m.exogenous_shocks))
            write(file, "shock_revtransforms", rev_transforms)
        end
    end
end

"""
```
write_forecast_block(file::JLD.JldFile, arr::Array, block_number::Int,
    block_inds::Range{Int64})
```

Writes `arr` to the subarray of `file` indicated by `block_inds`.
"""
function write_forecast_block(file::JLD.JldFile, arr::Array,
                              block_inds::Range{Int64})
    pfile = file.plain
    dataset = HDF5.d_open(pfile, "arr")
    ndims = length(size(dataset))
    dataset[block_inds, fill(Colon(), ndims-1)...] = arr
end

"""
```
read_forecast_metadata(file::JLD.JldFile)
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
function read_forecast_metadata(file::JLD.JldFile)
    metadata = Dict{Symbol, Any}()
    for field in setdiff(names(file), "arr")
        metadata[Symbol(field)] = read(file, field)
    end
    return metadata
end

"""
```
read_forecast_output(file, class, product, var_name[, shock_name])
```

Read only the forecast output for a particular variable (e.g. for a particular
observable) and possibly a particular shock. Result should be a matrix of size
`ndraws` x `nperiods`.
"""
function read_forecast_output(file::JLD.JldFile, class::Symbol, product::Symbol, var_name::Symbol)
    # Get index corresponding to var_name
    class_long = get_class_longname(class)
    indices = read(file, "$(class_long)_indices")
    var_ind = indices[var_name]

    pfile = file.plain
    filename = pfile.filename
    dataset = HDF5.d_open(pfile, "arr")
    ndims = length(size(dataset))

    # Trends are ndraws x nvars
    if product == :trend
        if ndims == 1 # one draw
            arr = h5read(filename, "arr", (var_ind,))
            arr = reshape(arr, (1, 1))
        elseif ndims == 2 # many draws
            arr = h5read(filename, "arr", (Colon(), var_ind))
        end

    # Other products are ndraws x nvars x nperiods
    elseif product in [:hist, :hist4q, :forecast, :bddforecast, :forecastut, :bddforecastut,
                       :forecast4q, :bddforecast4q, :dettrend]
        inds_to_read = if ndims == 2 # one draw
            arr = h5read(filename, "arr", (var_ind, Colon()))
        elseif ndims == 3 # many draws
            arr = h5read(filename, "arr", (Colon(), var_ind, Colon()))
            arr = squeeze(arr, 2)
        end
    else
        error("Invalid product: $product for this method")
    end

    return arr
end
function read_forecast_output(file::JLD.JldFile, class::Symbol, product::Symbol, var_name::Symbol,
                              shock_name::Symbol)
    # Get indices corresponding to var_name and shock_name
    class_long = get_class_longname(class)
    indices = read(file, "$(class_long)_indices")
    var_ind = indices[var_name]
    shock_indices = read(file, "shock_indices")
    shock_ind = shock_indices[shock_name]

    pfile = file.plain
    filename = pfile.filename
    dataset = HDF5.d_open(pfile, "arr")
    ndims = length(size(dataset))

    if ndims == 3 # one draw
        arr = h5read(filename, "arr", (var_ind, Colon(), shock_ind))
        arr = squeeze(arr, (1, 3))
        arr = reshape(arr, (1, length(arr)))
    elseif ndims == 4 # many draws
        arr = h5read(filename, "arr", (Colon(), var_ind, Colon(), shock_ind))
        arr = squeeze(arr, (2, 4))
    end

    return arr
end