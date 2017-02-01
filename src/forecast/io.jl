"""
```
get_forecast_input_file(m, input_type)
```

Compute the appropriate forecast input filenames for model `m` and
forecast input type `input_type`.

The default input files for each `input_type` can be overriden by adding entries
to the `Dict{Symbol, ASCIIString}` returned from
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
            error("Invalid input file override for input_type = $input_type: $override_file")
        end
    # If input_type = :subset or :block, also check for existence of overrides[:full]
    elseif input_type in [:subset, :block]
        return get_forecast_input_file(m, :full)
    end

    if input_type == :mode
        return rawpath(m,"estimate","paramsmode.h5")
    elseif input_type == :mean
        return workpath(m,"estimate","paramsmean.h5")
    elseif input_type == :init
        return ""
    elseif input_type in [:full, :subset, :block]
        return rawpath(m,"estimate","mhsave.h5")
    else
        throw(ArgumentError("Invalid input_type: $(input_type)"))
    end
end

"""
```
get_forecast_output_files(m, input_type, cond_type, output_vars;
                 [pathfcn = rawpath], [forecast_string = ""], [fileformat = :jld])
```

Compute the appropriate output filenames for model `m`, forecast input
type `input_type`, and conditional type `cond_type`, for each output variable in
`output_vars`. Returns a dictionary of file names with one entry for each output_var.

### Arguments

- See `forecast_one` for descriptions of other non-keyword arguments.

### Optional Arguments

- `pathfcn::Function`: should be one of `rawpath`, `workpath`,
  `inpath`, `figurespath`, `logpath`. Defaults to `rawpath`.
- `forecast_string::AbstractString`: subset identifier for when `input_type=:subset`
- `fileformat::Symbol`: file extension, without a period. Defaults to
  `:jld`, though `:h5` is another common option.


### Notes
See `get_forecast_filename` for more information.
"""
function get_forecast_output_files(m::AbstractModel,
                     input_type::Symbol, cond_type::Symbol, output_vars::Vector{Symbol};
                     forecast_string = "", fileformat = :jld)

    dir = rawpath(m, "forecast", "")
    model_string = get_model_string(m)

    get_forecast_output_files(dir, input_type, cond_type, output_vars,
                              model_string = model_string,
                              forecast_string = forecast_string,
                              fileformat = fileformat)
end

function get_forecast_output_files{S<:AbstractString}(directory::S,
                     input_type::Symbol, cond_type::Symbol, output_vars::Vector{Symbol};
                     model_string = "", forecast_string = "", fileformat = :jld)

    output_files = Dict{Symbol, S}()
    for var in output_vars
        output_files[var] = get_forecast_filename(directory, input_type, cond_type, var,
                                                  model_string=  model_string,
                                                  forecast_string = forecast_string,
                                                  fileformat = fileformat)
    end

    output_files
end

function get_model_string(m)

    tmp          = workpath(m, "forecast", "foo")
    output_dir   = dirname(tmp)
    replace(basename(tmp), "foo_", "")
end

"""
```
get_forecast_filename(m, input_type, cond_type, output_var;
                 [pathfcn = rawpath], [forecast_string = ""], [fileformat = :jld])


get_forecast_filename(directory, input_type, cond_type, output_var;
                 [forecast_string = ""], [fileformat = :jld])
```

### Notes

- If `input_type == :subset`, then the `forecast_string` is also appended to the
filenames. If in this case `forecast_string` is empty, `get_forecast_filename` throws
an error.

- In the second method, `directory` should be something of the form

```
\"\$saveroot/m990/ss2/forecast/raw/\"
```

Note that a `pathfcn` is therefore not required.

### Usage

```julia
output_file = get_forecast_filename(m, :mode, :none, :forecastpseudo)
```

The result will be something like:

```julia
\"\$saveroot/m990/ss2/forecast/raw/forecastpseudo_cond=none_para=mode_vint=REF.jld\"
```

Another example:

```julia
output_file = get_forecast_filename(m, :mode, :none, :forecastpseudo, workpath)
```

The result will be something like:

```julia
\"\$saveroot/m990/ss2/forecast/work/forecastpseudo_cond=none_para=mode_vint=REF.jld\"
```
"""
function get_forecast_filename(m::AbstractModel,
                     input_type::Symbol, cond_type::Symbol, output_var::Symbol;
                     pathfcn::Function = rawpath, model_string = "",
                     forecast_string = "", fileformat = :jld)

    # first, we need to make sure we get all of the settings that have been printed to this filestring
    fake_filename = pathfcn(m, "forecast", "foo")
    dir           = dirname(fake_filename)
    modelstring   = replace(basename(fake_filename), "foo_", "")

    # now we can find the filestring we are looking for
    get_forecast_filename(dir, input_type, cond_type, output_var,
                          model_string = model_string, forecast_string = forecast_string,
                          fileformat = fileformat)
end


function get_forecast_filename{S<:AbstractString}(directory::S,
                     input_type::Symbol, cond_type::Symbol, output_var::Symbol;
                     model_string = "", forecast_string = "", fileformat = :jld)

    # cast the strings so they're all the same type
    model_string = S(model_string)
    forecast_string = S(forecast_string)

    input_type = input_type == :block ? :full : input_type

    # gather all of the file strings into an array
    additional_file_strings = [model_string]
    push!(additional_file_strings, "para=" * abbrev_symbol(input_type))
    if isempty(forecast_string)
        if input_type == :subset
            error("Must supply a nonempty forecast_string if input_type = subset")
        end
    else
        push!(additional_file_strings, "fcid=" * forecast_string)
    end
    push!(additional_file_strings, "cond=" * abbrev_symbol(cond_type))

    sorted_model_string = filestring(additional_file_strings)
    joinpath(directory, "$output_var$(sorted_model_string).$fileformat")
end

"""
```
write_forecast_outputs(m, output_vars, forecast_output_files, forecast_output; verbose = :low)
```

Writes the elements of `forecast_output` indexed by `output_vars` to file, given
`forecast_output_files`, which maps `output_vars` to file names.
"""
function write_forecast_outputs{S<:AbstractString}(m::AbstractModel, input_type::Symbol,
                                output_vars::Vector{Symbol},
                                forecast_output_files::Dict{Symbol,S},
                                forecast_output::Dict{Symbol, Array{Float64}};
                                block_number::Nullable{Int64} = Nullable{Int64}(),
                                block_inds::Range{Int64} = 1:0,
                                verbose::Symbol = :low)

    for var in output_vars
        filepath = forecast_output_files[var]
        if isnull(block_number) || get(block_number) == 1
            jldopen(filepath, "w") do file
                write_forecast_metadata(m, file, var)
                write(file, "dims", get_forecast_output_dims(m, input_type, var))
            end
        end

        jldopen(filepath, "r+") do file
            if isnull(block_number)
                write(file, "arr", forecast_output[var])
            else
                write_forecast_block(file, forecast_output[var], get(block_number), block_inds)
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

- `date_indices::Dict{Date, Int}`: saved for all forecast outputs
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

        date_indices = [d::Date => i::Int for (i, d) in enumerate(dates)]
        write(file, "date_indices", date_indices)
    end

    # Write state names
    if class == :state
        state_indices = merge(m.endogenous_states, m.endogenous_states_augmented)
        @assert length(state_indices) == n_states_augmented(m) # assert no duplicate keys
        write(file, "state_indices", state_indices)
    end

    # Write observable names and transforms
    if class == :obs
        write(file, "observable_indices", m.observables)
        rev_transforms = if prod != :irf
            Dict{Symbol,Symbol}([x => symbol(m.observable_mappings[x].rev_transform) for x in keys(m.observables)])
        else
            Dict{Symbol,Symbol}([x => symbol("DSGE.identity") for x in keys(m.observables)])
        end
        write(file, "observable_revtransforms", rev_transforms)
    end

    # Write pseudo-observable names and transforms
    if class == :pseudo
        pseudo, pseudo_mapping = pseudo_measurement(m)
        write(file, "pseudoobservable_indices", pseudo_mapping.inds)
        rev_transforms = if prod != :irf
            Dict{Symbol,Symbol}([x => symbol(pseudo[x].rev_transform) for x in keys(pseudo)])
        else
            Dict{Symbol,Symbol}([x => symbol("DSGE.identity") for x in keys(pseudo)])
        end
        write(file, "pseudoobservable_revtransforms", rev_transforms)
    end

    # Write shock names and transforms
    if class == :shock || prod in [:shockdec, :irf]
        write(file, "shock_indices", m.exogenous_shocks)
        if class == :shock
            rev_transforms = Dict{Symbol,Symbol}([x => symbol("DSGE.identity") for x in keys(m.exogenous_shocks)])
        end
        write(file, "shock_revtransforms", rev_transforms)
    end
end

"""
```
write_forecast_block(file::JLD.JldFile, arr::Array, block_number::Int,
    block_inds::Range{Int64})
```

Writes `arr` as \"arr\$(block_number)\" and `block_inds` as
\"draws\$(block_number)\" to `file`, and updates `nblocks` and `ndraws`.
"""
function write_forecast_block(file::JLD.JldFile, arr::Array,
                              block_number::Int, block_inds::Range{Int64})

    # If datasets already exist, delete them
    exists(file, "arr$(block_number)")   ? delete!(file, "arr$(block_number)")   : nothing
    exists(file, "draws$(block_number)") ? delete!(file, "draws$(block_number)") : nothing

    write(file, "arr$(block_number)", arr)
    write(file, "draws$(block_number)", block_inds)

    # Update nblocks if necessary
    if exists(file, "nblocks")
        nblocks = read(file, "nblocks")
        delete!(file, "nblocks")
        write(file, "nblocks", block_number)
    else
        write(file, "nblocks", 1)
    end
end

"""
```
read_forecast_output(file::JLD.JldFile)
```

Returns the metadata dictionary and forecast output array in `file`.
"""
function read_forecast_output(file::JLD.JldFile)
    metadata = read_forecast_metadata(file)
    arr = if exists(file, "nblocks")
        read_forecast_blocks(file)
    else
        read(file, "arr")
    end

    return metadata, arr
end

"""
```
read_forecast_metadata(file::JLD.JldFile)
```

Read metadata from forecast output files. This includes dictionaries mapping dates, as well as state, observable,
pseudo-observable, and shock names, to their respective indices in the saved
forecast output array. The saved dictionaries include:

- `date_indices::Dict{Date, Int}`: saved for all forecast outputs
- `state_names::Dict{Symbol, Int}`: saved for `var in [:histstates, :forecaststates, :shockdecstates]`
- `observable_names::Dict{Symbol, Int}`: saved for `var in [:forecastobs, :shockdecobs]`
- `observable_revtransforms::Dict{Symbol, Symbol}`: saved identifiers for reverse transforms used for observables
- `pseudoobservable_names::Dict{Symbol, Int}`: saved for `var in [:histpseudo, :forecastpseudo, :shockdecpseudo]`
- `pseudoobservable_revtransforms::Dict{Symbol, Symbol}`: saved identifiers for reverse transforms used for pseudoobservables
- `shock_names::Dict{Symbol, Int}`: saved for `var in [:histshocks, :forecastshocks, :shockdecstates, :shockdecobs, :shockdecpseudo]`
"""
function read_forecast_metadata(file::JLD.JldFile)
    metadata = Dict{Symbol, Any}()
    for field in names(file)
        metadata[symbol(field)] = read(file, field)
    end
    return metadata
end

"""
```
read_forecast_blocks(file::JLD.JldFile)
```

Reads a JLD file written by `write_forecast_outputs` calling
`write_forecast_block`. Returns an `Array`.
"""
function read_forecast_blocks(file::JLD.JldFile)
    @assert exists(file, "nblocks") && exists(file, "dims")
    nblocks = read(file, "nblocks")
    dims    = read(file, "dims")
    ndims   = length(dims)

    arr = zeros(dims...)
    current_draw = 0
    for block = 1:nblocks
        block_draws  = read(file, "draws$block")
        ndraws_block = length(block_draws)
        block_draws_thin = (current_draw+1):(current_draw+ndraws_block)

        arr[block_draws_thin, fill(Colon(), ndims-1)...] = read(file, "arr$block")

        current_draw += ndraws_block
    end
    return arr
end