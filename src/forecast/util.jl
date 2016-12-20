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
reset_procs(m, procs, input_type = Nullable{Symbol}())
```

Reset `procs` to `[myid()]` if either:

- `input_type` is non-null and one of `[:init, :mode, :mean]`
- `use_parallel_workers(m) = false`
"""
function reset_procs(m::AbstractModel, procs::Vector{Int}, input_type::Nullable{Symbol} = Nullable{Symbol}())
    if procs != [myid()]
        if !isnull(input_type) && get(input_type) in [:init, :mode, :mean]
            warn("Processes $procs passed in, but input_type = $input_type. Using only process $(myid()) instead.")
            procs = [myid()]
        elseif !use_parallel_workers(m)
            warn("Processes $procs passed in, but use_parallel_workers(m) = false. Using only process $(myid()) instead.")
            procs = [myid()]
        end
    end
    return procs
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
get_input_file(m, input_type)
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
function get_input_file(m, input_type)
    overrides = forecast_input_file_overrides(m)
    if haskey(overrides, input_type)
        override_file = overrides[input_type]
        if ispath(override_file)
            return override_file
        else
            error("Invalid input file override for input_type = $input_type: $override_file")
        end
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
get_output_files(m, base, input_type, output_vars, cond_type;
                 [pathfcn = rawpath], [subset_string = ""], [fileformat = :jld])
```

Compute the appropriate output filenames for model `m`, forecast input
type `input_type`, and conditional type `cond_type`, for each output variable in
`output_vars`. Returns a dictionary of file names with one entry for each output_var.

### Arguments

- `base::AbstractString`: the output subdirectory corresponding to
  this general stage of the DSGE procedure. Examples could be `estimate` or `forecast`.
- See `forecast_one` for descriptions of other non-keyword arguments.

### Optional Arguments

- `pathfcn::Function`: should be one of `rawpath`, `workpath`,
  `inpath`, `figurespath`, `logpath`. Defaults to `rawpath`.
- `subset_string::AbstractString`: subset identifier for when `input_type=:subset`
- `fileformat::Symbol`: file extension, without a period. Defaults to
  `:jld`, though `:h5` is another common option.


### Notes
- If `input_type == :subset`, then the `subset_string` is also appended to the
filenames. If in this case `subset_string` is empty, `get_output_files` throws
an error.

- `base` can be any string, but is likely \"forecast\". An example use case is given below:

```julia
output_files = get_output_files(m, \"forecast\", :mode, [:forecastpseudo], :none)
```

The entry corresponding to the `:forecastpseudo` key will look something like:

```julia
\"$saveroot/m990/ss2/forecast/raw/forecastpseudo_cond=none_para=mode_vint=REF.jld\"
```

Another example:

```julia
output_files = get_output_files(m, \"forecast\", :mode, [:forecastpseudo], :none, workpath)
```

The entry corresponding to the `:forecastpseudo` key will look something like:

```julia
\"$saveroot/m990/ss2/forecast/work/forecastpseudo_cond=none_para=mode_vint=REF.jld\"
```
"""
function get_output_files{S<:AbstractString}(m::AbstractModel, base::S,
                     input_type::Symbol, output_vars::Vector{Symbol}, cond_type::Symbol;
                     pathfcn::Function = rawpath, subset_string::S = "",
                     fileformat = :jld)
    additional_file_strings = ASCIIString[]
    push!(additional_file_strings, "para=" * abbrev_symbol(input_type))
    if input_type == :subset
        if isempty(subset_string)
            error("Must supply a nonempty subset_string if input_type = subset")
        else
            push!(additional_file_strings, "sub=" * subset_string)
        end
    end
    push!(additional_file_strings, "cond=" * abbrev_symbol(cond_type))

    return convert(Dict{Symbol, ASCIIString},
                   [symbol(x) => pathfcn(m, base, "$x.$fileformat", additional_file_strings) for x in output_vars])
end

typealias DVector{T, A} DArray{T, 1, A}
typealias DMatrix{T, A} DArray{T, 2, A}

"""
```
dinit(T, dims)

dinit(T, dims...)
```

Initialize a `DArray{T, N, Array{T, N}}` (where `N = length(dims)`), distributed
over only the current process `myid()`.
"""
dinit(T, dims::Tuple) = DArray(I -> Array(T, map(length, I)), dims, [myid()])
dinit(T, dims...)     = dinit(T, dims)

"""
```
@time_verbose ex
```

A macro that calls `@time ex` if `VERBOSITY[verbose] >= VERBOSITY[:high]`, else
just calls `ex`.
"""
macro time_verbose(ex)
    quote
        local val = if VERBOSITY[verbose] >= VERBOSITY[:high]
            $(@time esc(ex))
        else
            $(esc(ex))
        end
        val
    end
end

"""
```
transplant_history(history, last_hist_period)
```

Remove the smoothed states, shocks, or pseudo-observables corresponding to
conditional data periods. This is necessary because when we forecast with
conditional data, we smooth beyond the last historical period.
"""
function transplant_history{T<:AbstractFloat}(history::DArray{T, 3},
    last_hist_period::Int,)

    return convert(DArray, history[1:end, 1:end, 1:last_hist_period])
end

"""
```
transplant_forecast(history, forecast, last_hist_period)
```

Transplant the smoothed states, shocks, or pseudo-observables corresponding to
conditional data periods from the history to the forecast.
"""
function transplant_forecast{T<:AbstractFloat}(history::DArray{T, 3},
    forecast::DArray{T, 3}, last_hist_period::Int)

    procs  = collect(history.pids)
    nprocs = length(procs)

    (ndraws, nvars, nfcastperiods) = size(forecast)
    ncondperiods = size(history, 3) - last_hist_period
    nperiods     = ncondperiods + nfcastperiods
    cond_range   = (last_hist_period + 1):(last_hist_period + ncondperiods)

    return DArray((ndraws, nvars, nperiods), procs, [nprocs, 1, 1]) do I
        draw_inds    = first(I)
        cond_draws   = convert(Array, history[draw_inds, :, cond_range])
        fcast_draws  = convert(Array, forecast[draw_inds, :, :])
        return cat(3, cond_draws, fcast_draws)
    end
end

"""
```
transplant_forecast_observables(histstates, forecastobs, systems, last_hist_period)
```

Transplant the observables implied by `histstates` corresponding to
conditional data periods from the history to the forecast.

This exists as a separate function from `transplant_forecast` because we don't
usually map the smoothed historical states into observables, since they would
just give us back the data. However, in the conditional data periods, we only
have data for a subset of observables, so we need to get the remaining
observables by mapping the smoothed states.
"""
function transplant_forecast_observables{T<:AbstractFloat}(histstates::DArray{T, 3},
    forecastobs::DArray{T, 3}, systems::DVector{System{T}}, last_hist_period::Int)

    procs  = collect(histstates.pids)
    nprocs = length(procs)

    (ndraws, nvars, nfcastperiods) = size(forecastobs)
    ncondperiods = size(histstates, 3) - last_hist_period
    nperiods     = ncondperiods + nfcastperiods
    cond_range   = (last_hist_period + 1):(last_hist_period + ncondperiods)

    return DArray((ndraws, nvars, nperiods), procs, [nprocs, 1, 1]) do I
        draw_inds = first(I)
        ndraws_local = length(draw_inds)

        # Map conditional period smoothed states to observables
        cond_draws = zeros(ndraws_local, nvars, ncondperiods)
        for i in draw_inds
            i_local = mod(i-1, ndraws_local) + 1
            states_i = convert(Array, slice(histstates, i, :, cond_range))
            cond_draws[i_local, :, :] = systems[i][:ZZ]*states_i .+ systems[i][:DD]
        end

        fcast_draws = convert(Array, forecastobs[draw_inds, :, :])
        return cat(3, cond_draws, fcast_draws)
    end
end

"""
```
write_forecast_outputs(m, output_vars, forecast_output_files, forecast_output; verbose = :low)
```

Writes the elements of `forecast_output` indexed by `output_vars` to file, given
`forecast_output_files`, which maps `output_vars` to file names. Calls
`write_darray`.
"""
function write_forecast_outputs(m::AbstractModel, output_vars::Vector{Symbol},
                                forecast_output_files::Dict{Symbol,ASCIIString},
                                forecast_output::Dict{Symbol, DArray{Float64}};
                                verbose::Symbol = :low)

    for var in output_vars
        filepath = forecast_output_files[var]
        jldopen(filepath, "w") do file
            write_forecast_metadata(m, file, var)
        end
        write_darray(filepath, forecast_output[var])

        if VERBOSITY[verbose] >= VERBOSITY[:high]
            println(" * Wrote $(basename(filepath))")
        end
    end
end

"""
```
prepare_forecast_inputs!(m, input_type, cond_type, output_vars;
    df = DataFrame(), systems = dinit(System, 0), kals = dinit(Kalman{S}, 0),
    subset_inds = Vector{Int}(), verbose = :none, procs = [myid()])
```

Check that the provided inputs `df`, `systems`, `states`, `subset_inds`, and
`procs` are well-formed with respect to the provided `input_type`, `cond_type`,
and `output_vars`. If an input is not provided to the function, it is loaded
using the appropriate getter function.

### Inputs

- `m::AbstractModel`: model object
- `input_type::Symbol`: See documentation for `forecast_all`. Defaults to
  `:mode`
- `cond_type::Symbol`: See documentation for `forecast_all`. Defaults to `:none`
- `subset_inds::Vector{Int}`: indices specifying the draws we want to use. See
  `forecast_one` for more detail
- `output_vars::Vector{Symbol}`: vector of desired output variables. See
  `forecast_one`

### Keyword Arguments

- `df::DataFrame`: historical data. If `cond_type in [:semi, :full]`, then the
   final row of `df` should be the period containing conditional data. If not
   provided, then `df` will be loaded using `load_data` with the appropriate
   `cond_type`
- `systems::DVector{System{Float64}}`: vector of `n_sim_forecast` many `System`
  objects, one for each draw. If not provided, will be loaded using
  `prepare_systems`
- `kals::DVector{Kalman{Float64}}`: vector of `n_sim_forecast` many `Kalman`
  objects. If not provided, will be loaded using `filter_all`
- `subset_inds::Vector{Int}`: indices specifying the draws we want to use. If
  `input_type` is not `subset`, `subset_inds` will be ignored
- `verbose::Symbol`: desired frequency of function progress messages printed to
  standard out. One of `:none`, `:low`, or `:high`
- `procs::Vector{Int}`: list of worker processes that have been
  previously added by the user. Defaults to `[myid()]`

### Outputs

- `df`
- `systems`
- `kals`
- `procs`
"""
function prepare_forecast_inputs!{S<:AbstractFloat}(m::AbstractModel{S},
    input_type::Symbol, cond_type::Symbol, output_vars::Vector{Symbol};
    df::DataFrame = DataFrame(),
    systems::DVector{System{S}} = dinit(System{S}, 0),
    kals::DVector{Kalman{S}} = dinit(Kalman{S}, 0),
    subset_inds::Vector{Int} = Vector{Int}(),
    verbose::Symbol = :none, procs::Vector{Int} = [myid()])

    # Set forecast_pseudoobservables properly
    for output in output_vars
        if contains(string(output), "pseudo")
            m <= Setting(:forecast_pseudoobservables, true)
            break
        end
    end

    # Load data if not provided
    if isempty(df)
        df = load_data(m; cond_type = cond_type, try_disk = true, verbose = :none)
    else
        @assert df[1, :date] == date_presample_start(m)
        @assert df[end, :date] == (cond_type == :none ? date_mainsample_end(m) : date_conditional_end(m))
    end

    # Compute systems and run Kalman filter if not provided
    if isempty(systems) || isempty(kals)
        procs = reset_procs(m, procs, Nullable(input_type))
        params = load_draws(m, input_type; subset_inds = subset_inds,
                            verbose = verbose, procs = procs)
        systems = prepare_systems(m, input_type, params; procs = procs)
        kals = filter_all(m, df, systems; cond_type = cond_type, procs = procs)
    else
        @assert length(systems) == length(kals)
        @assert DistributedArrays.procs(systems) == DistributedArrays.procs(kals) == procs
        if input_type == :subset
            @assert length(subset_inds) == length(systems)
        end
    end

    return df, systems, kals, procs
end

"""
```
write_darray(filepath::AbstractString, darr::DArray)
```

Write the contents of the `DArray` `darr` to the JLD file `filepath` without
converting `darr` back to an `Array` or copying local parts to the originator
process. `darr` can then be read back in as an `Array` using `read_darray`.

Each worker process `pid` writes its own local part and local indices to the
file as `arr\$pid` and `inds\$pid`. The dimensions `dims` of `darr` and the
vector of processes `pids` over which `darr` is distributed are also written to
`filepath` in order to facilitate reading back in.
"""
function write_darray{T<:AbstractFloat}(filepath::AbstractString, darr::DArray{T})
    function write_localpart(pid::Int)
        jldopen(filepath, "r+") do file
            write(file, "inds$pid", collect(localindexes(darr)))
            write(file, "arr$pid", localpart(darr))
        end
    end

    mode = isfile(filepath) ? "r+" : "w"
    jldopen(filepath, mode) do file
        write(file, "dims", darr.dims)
        write(file, "pids", collect(darr.pids))
    end

    for pid in darr.pids
        remotecall_wait(pid, write_localpart, pid)
        sleep(0.001) # Need this, or else HDF5 complains that the new object already exists
    end
end

"""
```
read_darray(file::JldFile)
```

Read the `DArray` saved to `file` by `write_darray`. Returns an `Array`.
"""
function read_darray(file::JLD.JldFile)
    dims = read(file, "dims")
    pids = read(file, "pids")

    out = zeros(dims...)
    for pid in pids
        inds = read(file, "inds$pid")
        out[inds...] = read(file, "arr$pid")
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
"""
function write_forecast_metadata(m::AbstractModel, file::JLD.JldFile, var::Symbol)

    # Write date range
    dates = if contains(string(var), "hist")
        quarter_range(date_mainsample_start(m), date_mainsample_end(m))
    elseif contains(string(var), "forecast")
        quarter_range(date_forecast_start(m), date_forecast_end(m))
    elseif contains(string(var), "shockdec")
        quarter_range(date_shockdec_start(m), date_shockdec_end(m))
    elseif contains(string(var), "trend") # trend and dettrend
        quarter_range(date_mainsample_start(m), date_forecast_end(m))
    elseif contains(string(var), "irf")
        NaN
    end

    if typeof(dates) == Array{Date,1}
        date_indices = [d::Date => i::Int for (i, d) in enumerate(dates)]
        write(file, "date_indices", date_indices)
    end

    # Write state names
    if contains(string(var), "states")
        state_indices = merge(m.endogenous_states, m.endogenous_states_augmented)
        @assert length(state_indices) == n_states_augmented(m) # assert no duplicate keys
        write(file, "state_indices", state_indices)
    end

    # Write observable names
    if contains(string(var), "obs")
        write(file, "observable_indices", m.observables)
        rev_transforms =
            Dict{Symbol,Symbol}([x => symbol(m.observable_mappings[x].rev_transform) for x in keys(m.observables)])
        write(file, "observable_revtransforms", rev_transforms)
    end

    # Write pseudo-observable names and transforms
    if contains(string(var), "pseudo")
        pseudo, pseudo_mapping = pseudo_measurement(m)
        write(file, "pseudoobservable_indices", pseudo_mapping.inds)
        rev_transforms = Dict{Symbol,Symbol}([x => symbol(pseudo[x].rev_transform) for x in keys(pseudo)])
        write(file, "pseudoobservable_revtransforms", rev_transforms)
    end

    # Write shock names
    if contains(string(var), "shocks") || contains(string(var), "shockdec") || contains(string(var), "irf")
        write(file, "shock_indices", m.exogenous_shocks)
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

    metadata
end


"""
```
compile_forecast_one(m, output_vars; verbose = :low, procs = [myid()])
```

Run `forecast_one` once to just-in-time compile it, before running a
full-distribution forecast. `compile_forecast_one` computes the minimum number
of draws necessary given the forecast thinning step size (`jstep`) and number of
processes (`nprocs`), i.e. `jstep * nprocs`, and runs a full forecast on those
draws.
"""
function compile_forecast_one(m, output_vars; verbose = :low, procs = [myid()])
    # Compute mininum number of draws for given jstep and procs
    jstep = get_setting(m, :forecast_jstep)
    nprocs = length(procs)
    min_draws = jstep * nprocs

    # Call forecast_one with input_type = :subset
    subset_inds = collect(1:min_draws)
    forecast_outputs = forecast_one(m, :subset, :none, output_vars;
                                    subset_inds = subset_inds, subset_string = "compile",
                                    verbose = verbose, procs = procs)

    # Delete output files
    output_files = get_output_files(m, "forecast", :subset, output_vars, :none, subset_string = "compile")
    map(rm, collect(values(output_files)))
end

"""
```
add_requisite_output_vars(m::AbstractModel, output_vars::Vector{Symbol})
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

    shockdec_outputs = Base.filter(output -> contains(string(output), "shockdec"), output_vars)
    if !isempty(shockdec_outputs)
        classes = [get_class(output) for output in shockdec_outputs]
        dettrend_vars = [symbol("dettrend$c") for c in classes]
        trend_vars = [symbol("trend$c") for c in classes]
        output_vars = unique(vcat(output_vars, dettrend_vars, trend_vars))
    end

    return output_vars
end