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

Returns a pair of `Vector{Range{Int64}}`s, `block_inds` and `block_inds_thin`,
each of length equal to the number of forecast blocks. `block_inds[i]` is the
range of indices for block `i` before thinning by `jstep` and
`block_inds_thin[i]` is the range after thinning.
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
    all_inds = start_ind:jstep:end_ind
    end_ind_thin = convert(Int64, floor(end_ind / jstep))

    # Make sure block_size is a multiple of jstep
    block_size = forecast_block_size(m)
    if block_size % jstep != 0
        error("forecast_block_size(m) must be a multiple of jstep = $jstep")
    end
    block_size_thin = convert(Int64, block_size / jstep)
    nblocks = convert(Int64, ceil(ndraws / block_size))

    # Fill in draw indices for each block
    block_inds      = Vector{Range{Int64}}(nblocks)
    block_inds_thin = Vector{Range{Int64}}(nblocks)
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

    return block_inds, block_inds_thin
end


"""
```
add_requisite_output_vars(output_vars::Vector{Symbol})
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
"""
function add_requisite_output_vars(output_vars::Vector{Symbol})

    # Add :bddforecast<class> if :forecast<class> is in output_vars
    forecast_outputs = Base.filter(output -> get_product(output) in [:forecast, :forecast4q], output_vars)
    if !isempty(forecast_outputs)
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
transplant_history(history, last_hist_period)
```

Remove the smoothed states, shocks, or pseudo-observables corresponding to
conditional data periods. This is necessary because when we forecast with
conditional data, we smooth beyond the last historical period.
"""
function transplant_history{T<:AbstractFloat}(history::Matrix{T},
    last_hist_period::Int)

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

    return hcat(condobs, forecastobs)
end

"""
```
standardize_shocks(shocks, QQ)
```

Normalize shocks by their standard deviations. Shocks with zero standard
deviation will be set to zero.
"""
function standardize_shocks{T<:AbstractFloat}(shocks::Matrix{T}, QQ::Matrix{T})
    stdshocks = shocks ./ sqrt(diag(QQ))

    zeroed_shocks = find(diag(QQ) .== 0)
    stdshocks[zeroed_shocks, :] = 0

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
function assemble_block_outputs(dicts::Vector{Dict{Symbol, Array{Float64}}})
    out = Dict{Symbol, Array{Float64}}()
    if !isempty(dicts)
        for var in keys(dicts[1])
            outputs  = map(dict -> reshape(dict[var], (1, size(dict[var])...)), dicts)
            out[var] = cat(1, outputs...)
        end
    end
    return out
end

"""
```
get_forecast_output_dims(m, input_type, output_var)
```

Returns the dimension of the forecast output specified by `input_type` and
`output_var`.
"""
function get_forecast_output_dims(m::AbstractModel, input_type::Symbol, output_var::Symbol;
                                  subset_inds::Range{Int64} = 1:0)
    prod  = get_product(output_var)
    class = get_class(output_var)

    ndraws = if input_type in [:mode, :mean, :init]
        1
    elseif input_type in [:full, :subset]
        _, block_inds_thin = forecast_block_inds(m, input_type; subset_inds = subset_inds)
        sum(map(length, block_inds_thin))
    end

    nvars = if class == :state
        n_states_augmented(m)
    elseif class == :obs
        n_observables(m)
    elseif class == :pseudo
        n_pseudoobservables(m)
    elseif class in [:shock, :stdshock]
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
        return (ndraws, nvars)
    elseif prod in [:hist, :forecast, :bddforecast, :dettrend]
        return (ndraws, nvars, nperiods)
    elseif prod in [:shockdec, :irf]
        nshocks = n_shocks_exogenous(m)
        return (ndraws, nvars, nperiods, nshocks)
    end
end
