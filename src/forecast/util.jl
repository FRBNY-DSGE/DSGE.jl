"""
```
get_jstep(m, n_sim)
```

Retrieve `forecast_jstep` setting (thinning step size for forecast
step) from `m.settings`. If `n_sim ==  1`, returns 1.
"""
function get_jstep(m::AbstractDSGEModel, n_sim::Int)
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
function n_forecast_draws(m::AbstractDSGEModel, input_type::Symbol)
    if input_type in [:mean, :mode, :init]
        return 1
    elseif input_type in [:full, :subset]
        input_file = get_forecast_input_file(m, input_type)
        if input_file[end-1:end] == "h5"
            draws = h5open(input_file, "r") do file
                if get_setting(m, :sampling_method) == :MH
                    dataset = isdefined(HDF5, :open_object) ? HDF5.open_object(file, "mhparams") : HDF5.o_open(file, "mhparams")
                elseif get_setting(m, :sampling_method) == :SMC
                    dataset = isdefined(HDF5, :open_object) ? HDF5.open_object(file, "smcparams") : HDF5.o_open(file, "smcparams")
                else
                    throw("Invalid :sampling_method setting specification.")
                end
                size(dataset)[1]
            end
        elseif input_file[end-3:end] == "jld2"
            draws = jldopen(input_file, "r") do file
                dataset = file["cloud"].particles
                size(dataset)[1]
            end
        end
        return draws
    elseif input_type == :prior || input_type == :mode_draw_shocks
        return 5000
    else
        throw(ArgumentError("Invalid input_type: $(input_type)"))
    end
end

"""
```
forecast_block_inds(m, input_type; subset_inds = 1:0)
```

Returns a pair of `Vector{AbstractRange{Int64}}`s, `block_inds` and `block_inds_thin`,
each of length equal to the number of forecast blocks. `block_inds[i]` is the
range of indices for block `i` before thinning by `jstep` and
`block_inds_thin[i]` is the range after thinning.
"""
function forecast_block_inds(m::AbstractDSGEModel, input_type::Symbol; subset_inds::AbstractRange{Int64} = 1:0, params::Vector{Float64} = Vector{Float64}(undef, 0))

    if !isempty(params)
        ndraws    = 1000
        jstep     = get_jstep(m, ndraws)
        start_ind = 1
        end_ind   = ndraws
    else
        if input_type == :full || input_type == :prior || input_type == :init_draw_shocks || input_type == :mode_draw_shocks
            set_ndraw = haskey(m.settings, :forecast_ndraws) ? get_setting(m, :forecast_ndraws) : 0
            if set_ndraw > 0
                ndraws = set_ndraw
            else
                ndraws = n_forecast_draws(m, :full)
            end
            jstep     = get_jstep(m, ndraws)
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
    end
    all_inds = start_ind:jstep:end_ind
    end_ind_thin = convert(Int64, floor(end_ind / jstep))

    # Make sure block_size is a multiple of jstep
    block_size = forecast_block_size(m)
    if end_ind_thin < block_size
        if input_type != :subset
            error("Number of draw ($(end_ind_thin)) divided by jstep ($(jstep)) is less than forecast block size.")
        else
            error("Last index of subset of draws ($(end_ind_thin)) divided by jstep ($(jstep)) is less than forecast block size.")
        end
    end
    if block_size % jstep != 0
        error("forecast_block_size(m) must be a multiple of jstep = $jstep")
    end
    block_size_thin = convert(Int64, block_size / jstep)
    nblocks = convert(Int64, ceil(ndraws / block_size))

    # Fill in draw indices for each block
    block_inds      = Vector{AbstractRange{Int64}}(undef, nblocks)
    block_inds_thin = Vector{AbstractRange{Int64}}(undef, nblocks)
    current_draw      = start_ind - 1
    current_draw_thin = start_ind - 1
    for i = 1:(nblocks-1)
        block_inds[i]      = (current_draw+jstep):jstep:(current_draw+block_size)
        block_inds_thin[i] = (current_draw_thin+1):(current_draw_thin+block_size_thin)
        current_draw      += block_size
        current_draw_thin += block_size_thin
    end
    block_inds[end]      = (current_draw+jstep):jstep:end_ind
    block_inds_thin[end] = (current_draw_thin+1):end_ind_thin

    # Check that block indices are all nonempty
    empty_blocks      = Int64[]
    empty_blocks_thin = Int64[]
    for (i, block_ind) in enumerate(block_inds)
        if isempty(block_ind)
            push!(empty_blocks, i)
        end
    end
    for (i, block_ind) in enumerate(block_inds_thin)
        if isempty(block_ind)
            push!(empty_blocks_thin, i)
        end
    end

    if !isempty(empty_blocks)
        block_inds = block_inds[empty_blocks]
    end
    if !isempty(empty_blocks_thin)
        block_inds_thin = block_inds_thin[empty_blocks_thin]
    end

    return block_inds, block_inds_thin
end


"""
```
add_requisite_output_vars(output_vars::Vector{Symbol}; bdd_fcast = true)
```

Based on the given `output_vars`, this function determines which
additional `output_vars` must be computed and stored for future
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

The option `bdd_fcast` can allow us to avoid adding `bdd` forecast output vars
if these are not wanted by the user.
"""
function add_requisite_output_vars(output_vars::Vector{Symbol}; bdd_fcast::Bool = true)
    # Add :bddforecast<class> if :forecast<class> is in output_vars
    forecast_outputs = Base.filter(output -> get_product(output) in [:forecast, :forecastut, :forecast4q, :forecastlvl],
                                   output_vars)
    if !isempty(forecast_outputs) && bdd_fcast
        bdd_vars = [Symbol("bdd$(var)") for var in forecast_outputs]
        output_vars = unique(vcat(output_vars, bdd_vars))
    end

    # Add :trend<class> and :dettrend<class> if :shockdec<class> is in output_vars
    shockdec_outputs = Base.filter(output -> get_product(output) == :shockdec, output_vars)
    if !isempty(shockdec_outputs)
        classes = [get_class(output) for output in shockdec_outputs]
        dettrend_vars = [Symbol("dettrend$c") for c in classes]
        trend_vars = [Symbol("trend$c") for c in classes]
        output_vars = unique(vcat(output_vars, dettrend_vars, trend_vars))
    end

    return output_vars
end

"""
```
remove_meansbands_only_output_vars(output_vars)
```
"""
function remove_meansbands_only_output_vars(output_vars::Vector{Symbol})
    # All the <product>ut<class> and <product>4q<class> variables are computed
    # during compute_meansbands
    meansbands_only_products = [:histut, :hist4q, :forecastut, :forecast4q,
                                :bddforecastut, :bddforecast4q, :histlvl, :forecastlvl,
                                :bddforecastlvl]

    Base.filter(var -> !(get_product(var) in meansbands_only_products), output_vars)
end

"""
```
transplant_history(history, last_hist_period)
```

Remove the smoothed states, shocks, or pseudo-observables corresponding to
conditional data periods. This is necessary because when we forecast with
conditional data, we smooth beyond the last historical period.
"""
function transplant_history(history::Matrix{T}, last_hist_period::Int) where {T<:AbstractFloat}

    if isempty(history)
        return history
    else
        return history[:, 1:last_hist_period]
    end
end

"""
```
transplant_forecast(history, forecast, last_hist_period)
```

Transplant the smoothed states, shocks, or pseudo-observables corresponding to
conditional data periods from the history to the forecast.
"""
function transplant_forecast(history::Matrix{T}, forecast::Matrix{T},
                             last_hist_period::Int) where {T<:AbstractFloat}

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
function transplant_forecast_observables(histstates::Matrix{T}, forecastobs::Matrix{T},
                                         system::System{T},
                                         last_hist_period::Int) where {T<:AbstractFloat}

    nvars        = size(forecastobs, 1)
    ncondperiods = size(histstates, 2) - last_hist_period
    cond_range   = (last_hist_period + 1):(last_hist_period + ncondperiods)

    condstates   = histstates[:, cond_range]
    condobs      = system[:ZZ]*condstates .+ system[:DD]

    return hcat(condobs, forecastobs)
end

"""
```
standardize_shocks(shocks, QQ)
```

Normalize shocks by their standard deviations. Shocks with zero standard
deviation will be set to zero.
"""
function standardize_shocks(shocks::AbstractMatrix{T}, QQ::AbstractMatrix{T}) where {T<:AbstractFloat}
    stdshocks = shocks ./ sqrt.(diag(QQ))

    zeroed_shocks = findall(diag(QQ) .== 0)
    stdshocks[zeroed_shocks, :] .= 0

    return stdshocks
end

"""
```
standardize_shocks(shocks, QQs, start_date, end_date)
```

Normalize shocks by their standard deviations when there is regime switching.
Shocks with zero standard deviation will be set to zero.
"""
function standardize_shocks(shocks::Matrix{T}, QQs::Vector{Matrix{T}},
                            regime_inds::Vector{UnitRange{Int}}) where {T<:AbstractFloat}

    stdshocks = similar(shocks)
    for(QQ, inds) in zip(QQs, regime_inds)
        inshocks = @view shocks[:, inds]
        stdshocks[:, inds] = standardize_shocks(inshocks, QQ)
    end

    return stdshocks
end

"""
```
assemble_block_outputs(dicts)
```

Given a vector `dicts` of forecast output dictionaries, concatenate each output
along the draw dimension and return a new dictionary of the concatenated
outputs.
"""
function assemble_block_outputs(dicts::Vector{Dict{Symbol, Array{Float64}}}; show_failed_percent::Bool = false)
    out = Dict{Symbol, Array{Float64}}()
    if !isempty(dicts)
        _populate_empty_dictionaries!(dicts; show_failed_percent = show_failed_percent)
        for var in keys(dicts[1])
            outputs  = map(dict -> reshape(dict[var], (1, size(dict[var])...)), dicts)
            out[var] = cat(outputs..., dims=1)
        end
    end
    return out
end

"""
```
_populate_empty_dictionaries!(dicts)
```

helps `assemble_block_outputs` handle empty dictionaries by populating
them with NaNs.
"""
function _populate_empty_dictionaries!(dicts::Vector{Dict{Symbol, Array{Float64}}}; show_failed_percent::Bool = false)
    check = isempty.(dicts)
    if any(check)
        empty_dicts  = findall(check)
        if length(empty_dicts) == length(dicts)
            error("All dictionaries in dicts are empty. No valid forecast output was computed by forecast_one_draw")
        else
            ref_dict = dicts[findfirst(.!check)]
            for i in empty_dicts
                for (k, v) in ref_dict
                    dicts[i][k] = fill(NaN, size(v))
                end
            end
            if show_failed_percent
                nan_percent = round(100. * count(check) / length(dicts), digits = 2)
                println("The percentage of failed forecasts is $(nan_percent)%")
            end
        end
    end

    ct = 0
    for dict in dicts
        for (k, v) in dict
            if any(isnan.(v))
                ct += 1
                break
            end
        end
    end
    if show_failed_percent
        nan_percent = round(100. * ct / length(dicts), digits = 2)
        println("The percentage of failed forecasts is $(nan_percent)%")
    end

end

"""
```
get_forecast_output_dims(m, input_type, output_var)
```

Returns the dimension of the forecast output specified by `input_type` and
`output_var`.
"""
function get_forecast_output_dims(m::AbstractDSGEModel, input_type::Symbol, output_var::Symbol;
                                  subset_inds::AbstractRange{Int64} = 1:0)
    prod  = get_product(output_var)
    class = get_class(output_var)

    ndraws = if input_type in [:mode, :mean, :init]
        1
    elseif input_type in [:full, :subset, :prior, :init_draw_shocks, :mode_draw_shocks]
        _, block_inds_thin = forecast_block_inds(m, input_type; subset_inds = subset_inds)
        sum(map(length, block_inds_thin))
    end

    nvars = if class == :states
        n_states_augmented(m)
    elseif class == :obs
        n_observables(m)
    elseif class == :pseudo
        n_pseudo_observables(m)
    elseif class in [:shocks, :stdshocks]
        n_shocks_exogenous(m)
    end

    nperiods = if prod == :hist
        n_mainsample_periods(m)
    elseif prod in [:forecast, :bddforecast]
        forecast_horizons(m)
    elseif prod in [:shockdec, :dettrend]
        n_shockdec_periods(m)
    elseif prod == :irf
        impulse_response_horizons(m)
    end

    if prod == :trend
        regime_switching = haskey(m.settings, :regime_switching) ? get_setting(m, :regime_switching) : false
        if regime_switching
            if haskey(get_settings(m), :time_varying_trends) && get_setting(m, :time_varying_trends)
                start_date = get(get_setting(m, :shockdec_startdate))
                end_date = max(prev_quarter(date_forecast_start(m)), date_forecast_end(m))
                return (ndraws, nvars, DSGE.subtract_quarters(end_date, start_date)+1)
            else regime_switching
                return (ndraws, nvars, get_setting(m, :n_regimes))
            end
        else
            return (ndraws, nvars)
        end
    elseif prod in [:hist, :forecast, :bddforecast, :dettrend]
        return (ndraws, nvars, nperiods)
    elseif prod in [:shockdec, :irf]
        nshocks = n_shocks_exogenous(m)
        return (ndraws, nvars, nperiods, nshocks)
    end
end
