"""
```
compute_system(m)
```

Given the current model parameters, compute the state-space system
corresponding to model `m`. Returns a `System` object.
"""
function compute_system{T<:AbstractFloat}(m::AbstractModel{T})
    # Solve model
    TTT, RRR, CCC = solve(m)
    transition_equation = Transition(TTT, RRR, CCC)

    # Solve measurement equation
    shocks = n_anticipated_shocks(m) > 0
    measurement_equation = measurement(m, TTT, RRR, CCC; shocks = shocks)

    # Solve pseudo-measurement equation
    pseudo_measurement_equation = if method_exists(pseudo_measurement, (typeof(m),)) && forecast_pseudoobservables(m)
        _, pseudo_mapping = pseudo_measurement(m)
        Nullable(pseudo_mapping)
    else
        Nullable{PseudoObservableMapping{T}}()
    end

    return System(transition_equation, measurement_equation, pseudo_measurement_equation)
end

"""
```
get_jstep(m, n_sim)
```

Retrieve `forecast_jstep` setting (thinning step size for forecast
step) from `m.settings`. If `n_sim ==  1`, returns 1.
"""
function get_jstep(m, n_sim)
    if n_sim == 1
        jstep = 1
    else
        jstep = get_setting(m, :forecast_jstep)
    end
end

"""
```
n_forecast_draws(m, input_type)
```

Returns the number of forecast draws in the file
`get_forecast_input_file(m, input_type)`.
"""
function n_forecast_draws(m::AbstractModel, input_type::Symbol)
    if input_type in [:mean, :mode, :init]
        return 1
    elseif input_type in [:full, :subset]
        input_file = get_forecast_input_file(m, input_type)
        draws = h5open(input_file, "r") do file
            dataset = HDF5.o_open(file, "mhparams")
            size(dataset)[1]
        end
        return draws
    else
        throw(ArgumentError("Invalid input_type: $(input_type)"))
    end
end

"""
```
forecast_block_inds(m, input_type; subset_inds = 1:0)
```

Returns a `Vector{Range{Int64}}` of length equal to the number of forecast
blocks, where `block_inds[i]` is the range of indices (after thinning by
`jstep`) for block `i`.
"""
function forecast_block_inds(m::AbstractModel, input_type::Symbol; subset_inds::Range{Int64} = 1:0)

    if input_type == :full
        ndraws = n_forecast_draws(m, :full)
        jstep = get_jstep(m, ndraws)
        start_ind = 1
        end_ind   = ndraws
    elseif input_type == :subset
        ndraws    = length(subset_inds)
        jstep     = get_jstep(m, ndraws)
        start_ind = first(subset_inds)
        end_ind   = last(subset_inds)
    else
        throw(ArgumentError("Cannot call forecast_block_inds with input_type = $input_type."))
    end
    all_inds  = start_ind:jstep:end_ind

    # Make sure block_size is a multiple of jstep
    block_size = forecast_block_size(m)
    if block_size % jstep != 0
        error("forecast_block_size(m) must be a multiple of jstep = $jstep")
    end
    nblocks = convert(Int64, ceil(ndraws / block_size))

    # Fill in draw indices for each block
    block_inds  = Vector{Range{Int64}}(nblocks)
    current_draw = start_ind - 1
    for i = 1:(nblocks-1)
        block_inds[i] = (current_draw+jstep-1):jstep:(current_draw+block_size)
        current_draw += block_size
    end
    block_inds[end] = (current_draw+1):end_ind

    return block_inds
end

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
transplant_history(history, last_hist_period)
```

Remove the smoothed states, shocks, or pseudo-observables corresponding to
conditional data periods. This is necessary because when we forecast with
conditional data, we smooth beyond the last historical period.
"""
function transplant_history{T<:AbstractFloat}(history::Matrix{T},
    last_hist_period::Int)

    return history[:, 1:last_hist_period]
end

"""
```
transplant_forecast(history, forecast, last_hist_period)
```

Transplant the smoothed states, shocks, or pseudo-observables corresponding to
conditional data periods from the history to the forecast.
"""
function transplant_forecast{T<:AbstractFloat}(history::Matrix{T},
    forecast::Matrix{T}, last_hist_period::Int)

    ncondperiods = size(history, 2) - last_hist_period
    cond_range   = (last_hist_period + 1):(last_hist_period + ncondperiods)
    condhist     = history[:, cond_range]

    return hcat(condhist, forecast)
end

"""
```
transplant_forecast_observables(histstates, forecastobs, system, last_hist_period)
```

Transplant the observables implied by `histstates` corresponding to
conditional data periods from the history to the forecast.

This exists as a separate function from `transplant_forecast` because we don't
usually map the smoothed historical states into observables, since they would
just give us back the data. However, in the conditional data periods, we only
have data for a subset of observables, so we need to get the remaining
observables by mapping the smoothed states.
"""
function transplant_forecast_observables{T<:AbstractFloat}(histstates::Matrix{T},
    forecastobs::Matrix{T}, system::System{T}, last_hist_period::Int)

    nvars        = size(forecastobs, 1)
    ncondperiods = size(histstates, 2) - last_hist_period
    cond_range   = (last_hist_period + 1):(last_hist_period + ncondperiods)

    condstates   = histstates[:, cond_range]
    condobs      = system[:ZZ]*condstates .+ system[:DD]

    return hcat(cond_draws, forecastobs)
end

"""
```
write_forecast_outputs(m, output_vars, forecast_output_files, forecast_output; verbose = :low)
```

Writes the elements of `forecast_output` indexed by `output_vars` to file, given
`forecast_output_files`, which maps `output_vars` to file names.
"""
function write_forecast_outputs{S<:AbstractString}(m::AbstractModel, output_vars::Vector{Symbol},
                                forecast_output_files::Dict{Symbol,S},
                                forecast_output::Dict{Symbol, Array{Float64}};
                                block_number::Nullable{Int64} = Nullable{Int64}(),
                                subset_inds::Range{Int64} = 1:0,
                                verbose::Symbol = :low)

    for var in output_vars
        filepath = forecast_output_files[var]
        if isnull(block_number) || get(block_number) == 1
            jldopen(filepath, "w") do file
                write_forecast_metadata(m, file, var)
            end
        end

        jldopen(filepath, "r+") do file
            if isnull(block_number)
                write(file, "arr", forecast_output[var])
            else
                write_forecast_block(file, forecast_output[var], get(block_number), subset_inds)
            end
        end

        if VERBOSITY[verbose] >= VERBOSITY[:high]
            println(" * Wrote $(basename(filepath))")
        end
    end
end

"""
```
assemble_block_outputs(dicts)
```

Given a vector `dicts` of forecast output dictionaries, concatenate each output
along the draw dimension and return a new dictionary of the concatenated
outputs.
"""
function assemble_block_outputs(dicts::Vector{Dict{Symbol, Array{Float64}}})
    out = Dict{Symbol, Array{Float64}}()
    if !isempty(dicts)
        for var in keys(dicts[1])
            outputs  = map(dict -> dict[var], dicts)
            out[var] = cat(1, outputs...)
        end
    end
    return out
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

    var = string(var)

    # Write date range
    if !contains(var, "irf")
        dates = if contains(var, "hist")
            quarter_range(date_mainsample_start(m), date_mainsample_end(m))
        elseif contains(var, "forecast")
            quarter_range(date_forecast_start(m), date_forecast_end(m))
        elseif contains(var, "shockdec") || contains(var, "trend") # trend and dettrend
            quarter_range(date_shockdec_start(m), date_shockdec_end(m))
        end

        date_indices = [d::Date => i::Int for (i, d) in enumerate(dates)]
        write(file, "date_indices", date_indices)
    end

    # Write state names
    if contains(var, "states")
        state_indices = merge(m.endogenous_states, m.endogenous_states_augmented)
        @assert length(state_indices) == n_states_augmented(m) # assert no duplicate keys
        write(file, "state_indices", state_indices)
    end

    # Write observable names and transforms
    if contains(var, "obs")
        write(file, "observable_indices", m.observables)
        rev_transforms = if !contains(var, "irf")
            Dict{Symbol,Symbol}([x => symbol(m.observable_mappings[x].rev_transform) for x in keys(m.observables)])
        else
            Dict{Symbol,Symbol}([x => symbol("DSGE.identity") for x in keys(m.observables)])
        end
        write(file, "observable_revtransforms", rev_transforms)
    end

    # Write pseudo-observable names and transforms
    if contains(var, "pseudo")
        pseudo, pseudo_mapping = pseudo_measurement(m)
        write(file, "pseudoobservable_indices", pseudo_mapping.inds)
        rev_transforms = if !contains(var, "irf")
            Dict{Symbol,Symbol}([x => symbol(pseudo[x].rev_transform) for x in keys(pseudo)])
        else
            Dict{Symbol,Symbol}([x => symbol("DSGE.identity") for x in keys(pseudo)])
        end
        write(file, "pseudoobservable_revtransforms", rev_transforms)
    end

    # Write shock names and transforms
    if contains(var, "shocks") || contains(var, "shockdec") || contains(var, "irf")
        write(file, "shock_indices", m.exogenous_shocks)
        if contains(var, "shocks")
            rev_transforms = Dict{Symbol,Symbol}([x => symbol("DSGE.identity") for x in keys(m.exogenous_shocks)])
        end
        write(file, "shock_revtransforms", rev_transforms)
    end
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
write_forecast_block(file::JLD.JldFile, arr::Array, block_number::Int,
    subset_inds::Range{Int64})
```

Writes `arr` as \"arr\$(block_number)\" and `subset_inds` as
\"draws\$(block_number)\" to `file`, and updates `nblocks` and `ndraws`.
"""
function write_forecast_block(file::JLD.JldFile, arr::Array,
                              block_number::Int, subset_inds::Range{Int64})

    write(file, "arr$(block_number)", arr)
    write(file, "draws$(block_number)", subset_inds)

    # Update nblocks and dims if necessary
    if exists(file, "nblocks")
        nblocks = read(file, "nblocks")
        delete!(file, "nblocks")
        write(file, "nblocks", nblocks + 1)
    else
        write(file, "nblocks", 1)
    end

    if exists(file, "dims")
        dims = read(file, "dims")
        dims[1] += length(subset_inds)
        delete!(file, "dims")
        write(file, "dims", dims)
    else
        write(file, "dims", collect(size(arr)))
    end
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
    for block = 1:nblocks
        block_draws = read(file, "draws$block")
        arr[block_draws, fill(:, ndims-1)...] = read(file, "arr$block")
    end
    return arr
end

"""
```
add_requisite_output_vars(output_vars::Vector{Symbol})
```

Based on the given output_vars, this function determines which
additional output_vars must be computed and stored for future
plotting.

Specifically, when plotting a shock decomposition, the trend and
deterministic trend series are required (the trend is subtracted from
the value of each shock, and the deterministic trend represents
deviations from steady-state that would realize even in the absence of
shocks). For example, if `output_vars` contains `shockdecobs`, the
variables `dettrendobs` and `trendobs` will be added to `output_vars`.

Note that this case is distinct from a case in which computing a
different product is required to compute the desired `output_var`. For
example, smoothed historical states (`histstates`) must be computed in
order to compute a shock decomposition for a state variable, but need
not be saved to produce plots later on. Therefore, `histstates` is not
added to `output_vars` when calling
`add_requisite_output_vars([shockdecstates])`.
"""
function add_requisite_output_vars(output_vars::Vector{Symbol})

    # Add :forecast<class>bdd if :forecast<class> is in output_vars
    forecast_outputs = Base.filter(output -> contains(string(output), "forecast") && !contains(string(output), "bdd"),
                                   output_vars)
    if !isempty(forecast_outputs)
        bdd_vars = [symbol("bdd$(var)") for var in forecast_outputs]
        output_vars = unique(vcat(output_vars, bdd_vars))
    end

    # Add :trend<class> and :dettrend<class> if :shockdec<class> is in output_vars
    shockdec_outputs = Base.filter(output -> contains(string(output), "shockdec"), output_vars)
    if !isempty(shockdec_outputs)
        classes = [get_class(output) for output in shockdec_outputs]
        dettrend_vars = [symbol("dettrend$c") for c in classes]
        trend_vars = [symbol("trend$c") for c in classes]
        output_vars = unique(vcat(output_vars, dettrend_vars, trend_vars))
    end

    return output_vars
end