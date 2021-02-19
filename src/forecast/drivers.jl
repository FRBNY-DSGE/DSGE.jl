"""
```
prepare_forecast_inputs!(m, input_type, cond_type, output_vars;
    df = DataFrame(), verbose = :none)
```

Add required outputs using `add_requisite_output_vars` and load data if
necessary.

### Inputs

- `m::AbstractDSGEModel`: model object
- `input_type::Symbol`: See `?forecast_one`.
- `cond_type::Symbol`: See `?forecast_one`.
- `output_vars::Vector{Symbol}`: vector of desired output variables. See
  `?forecast_one_draw`

### Keyword Arguments

- `df::DataFrame`: historical data. If `cond_type in [:semi, :full]`, then the
   final row of `df` should be the period containing conditional data. If not
   provided, then `df` will be loaded using `load_data` with the appropriate
   `cond_type`
- `only_filter::Bool`: do not run the smoother and only run the filter. This limits the number of
  output variables which can be calculated.
- `verbose::Symbol`: desired frequency of function progress messages printed to
  standard out. One of `:none`, `:low`, or `:high`

### Outputs

- `output_vars`
- `df`
"""
function prepare_forecast_inputs!(m::AbstractDSGEModel{S},
    input_type::Symbol, cond_type::Symbol, output_vars::Vector{Symbol};
    df::DataFrame = DataFrame(), subset_inds::AbstractRange{Int64} = 1:0,
    check_empty_columns::Bool = true, bdd_fcast::Bool = true, only_filter::Bool = false,
    verbose::Symbol = :none) where {S<:AbstractFloat}

    @assert cond_type in [:none, :semi, :full] "cond_type must be one of :none, :semi, or :full"

    if only_filter
        smooth_vars = [:histstates, :histpseudo, :histshocks, :histstdshocks, :forecastshocks, :forecaststdshocks,
                       :shockdecstates, :shockdecpseudo, :shockdecobs, :dettrendstates, :dettrendpseudo, :dettrendobs]
        @assert isempty(intersect(output_vars, smooth_vars)) "Some output variables requested cannot be use with keyword only_filter"
    end

    implied_horizons = subtract_quarters(date_forecast_end(m), date_forecast_start(m)) + 1
    if implied_horizons != get_setting(m, :forecast_horizons)
        @warn "Setting :forecast_horizons does not match the number of periods implied by Settings :date_forecast_end" *
            " and :date_forecast_start. Updating :forecast_horizons to match"
        m <= Setting(:forecast_horizons, implied_horizons, "Number of periods to forecast ahead")
    end

    # Compute everything that will be needed to plot original output_vars
    output_vars = add_requisite_output_vars(output_vars, bdd_fcast = bdd_fcast)
    if input_type == :prior
        output_vars = setdiff(output_vars, [:bddforecastobs])
    end

    # Get products and classes computed
    output_prods   = unique(map(get_product, output_vars))
    output_classes = unique(map(get_class,   output_vars))

    # Throw error if input_type = :subset but no subset_inds provided
    if input_type == :subset && isempty(subset_inds)
        error("Must supply nonempty subset_inds if input_type = :subset")
    end

    # Determine if we are only running IRFs. If so, we won't need to load data below
    irfs_only = all(prod -> prod == :irf, output_prods)

    # Load data if not provided, else check data well-formed
    if !irfs_only
        if isempty(df)
            data_verbose = verbose == :none ? :none : :low
            df = load_data(m; cond_type = cond_type, try_disk = true, check_empty_columns = check_empty_columns, verbose = data_verbose)
        else
            @assert df[1, :date] == date_presample_start(m)
            @assert df[end, :date] == (cond_type == :none ? date_mainsample_end(m) : date_conditional_end(m))
        end
    end

    # If regime switching, check settings are consistent
    regswitch = haskey(get_settings(m), :regime_switching) ? get_setting(m, :regime_switching) : false
    if regswitch
        @assert get_setting(m, :n_regimes) == length(get_setting(m, :regime_dates)) "The number" *
            " in Setting :n_regimes ($(string(get_setting(m, :n_regimes))))" *
            " does not match the length of Setting :regime_dates ($(string(length(get_setting(m, :regime_dates)))))."
        @assert get_setting(m, :regime_dates)[1] == date_presample_start(m) "The first regime" *
            " date ($(string(get_setting(m, :regime_dates)[1]))) must match the presample start date " *
            "($(string(date_presample_start(m))))."
    end

    # Print out settings of temporary policies
    if verbose == :high && regswitch
        @info "Regime-Switching and Temporary Policy Settings"
        println("n_regimes = $(get_setting(m, :n_regimes))")
        println("reg_forecast_start = $(get_setting(m, :reg_forecast_start))")
        println("reg_post_conditional_end = $(get_setting(m, :reg_post_conditional_end))")
        println("time-varying information set: ", haskey(get_settings(m), :tvis_information_set) ?
                get_setting(m, :tvis_information_set) : false)

        println("gensys2: " * string(haskey(get_settings(m), :gensys2) ? get_setting(m, :gensys2) : false))
        replace_eqcond = haskey(get_settings(m), :replace_eqcond) ? get_setting(m, :replace_eqcond) : false
        println("replace_eqcond: " * string(replace_eqcond))
        if replace_eqcond
            replace_regs = sort!(collect(keys(get_setting(m, :regime_eqcond_info))))
            println("regime_eqcond_info: ", OrderedDict(k => get_setting(m, :regime_eqcond_info)[k] for k in replace_regs))
        end

        println("uncertain_altpolicy: " *
                string(haskey(get_settings(m), :uncertain_altpolicy) ? get_setting(m, :uncertain_altpolicy) : false))
        println("uncertain_temp_altpol: " * string(haskey(get_settings(m), :uncertain_temp_altpol) ? get_setting(m, :uncertain_temp_altpol) : false))
        println("time-varying credibility: " * string(any(v -> !ismissing(v.weights), values(get_setting(m, :regime_eqcond_info)))))
        if haskey(get_settings(m), :uncertain_altpolicy) ? get_setting(m, :uncertain_altpolicy) : false
            println("Desired policy rule: " *
                    string(haskey(get_settings(m), :regime_eqcond_info) ?
                           get_setting(m, :regime_eqcond_info)[maximum(collect(keys(get_setting(m, :regime_eqcond_info))))] :
                           alternative_policy(m).key))
            println("Other policy rules: " * join([string(x.key) for x in get_setting(m, :alternative_policies)], ", "))
            if haskey(get_settings(m), :imperfect_awareness_weights) # does not work currently, need to update
                println("Credibility weights: ", get_setting(m, :imperfect_awareness_weights))
            else
                sorted_eqcond_info_regs = sort!(collect(keys(get_setting(m, :regime_eqcond_info))))
                println("Credbility weights: ", OrderedDict(reg => get_setting(m, :regime_eqcond_info)[reg].weights for reg in sorted_eqcond_info_regs))
            end
        end
        if haskey(get_settings(m), :temporary_altpol_length)
            println("temporary altpol length = ", get_setting(m, :temporary_altpol_length))
        end
        altkey = alternative_policy(m).key
        if altkey in [:smooth_ait_gdp_alt, :smooth_ait_gdp, :flexible_ait]
            println("skip_altpolicy_state_init: ", haskey(get_settings(m), :skip_altpolicy_state_init) ?
                get_setting(m, :skip_altpolicy_state_init) : false)
            println("add_initialize_pgap_ygap_pseudoobs: ", haskey(get_settings(m), :add_initialize_pgap_ygap_pseudoobs) ?
                get_setting(m, :add_initialize_pgap_ygap_pseudoobs) : false)
            println("φ_π = ", get_setting(m, Symbol(altkey, "_φ_π")))
            println("φ_y = ", get_setting(m, Symbol(altkey, "_φ_y")))
            println("ait_Thalf = ", get_setting(m, :ait_Thalf))
            println("gdp_Thalf = ", get_setting(m, :gdp_Thalf))
            println("ρ_smooth = ", get_setting(m, Symbol(altkey, "_ρ_smooth")))
            println("pgap type = ", get_setting(m, :pgap_type))
            println("ygap type = ", get_setting(m, :ygap_type))
            println("initial pgap = ", get_setting(m, :pgap_value))
            println("initial ygap = ", get_setting(m, :ygap_value))
        end

        # TODO: print out the forecast_zlb_value used by zero_rate
        println("") # add some space
    end

    return output_vars, df
end

"""
```
load_draws(m, input_type; subset_inds = 1:0, verbose = :low)

load_draws(m, input_type, block_inds; verbose = :low)
```

Load and return parameter draws from Metropolis-Hastings or SMC.

### Inputs

- `m::AbstractDSGEModel`: model object
- `input_type::Symbol`: one of the options for `input_type` described in the
  documentation for `forecast_one`
- `block_inds::AbstractRange{Int64}`: indices of the current block (already indexed by
  `jstep`) to be read in. Only used in second method

### Keyword Arguments

- `subset_inds::AbstractRange{Int64}`: indices specifying the subset of draws to be read
  in. Only used in first method
- `verbose::Symbol`: desired frequency of function progress messages printed to
  standard out. One of `:none`, `:low`, or `:high`. If `:low` or greater, prints
  location of input file.

### Outputs

- `params`: first method returns a single parameter draw of type
  `Vector{Float64}`. Second method returns a `Vector{Vector{Float64}}` of
  parameter draws for this block.
"""
function load_draws(m::AbstractDSGEModel, input_type::Symbol;
                    subset_inds::AbstractRange{Int64} = 1:0,
                    verbose::Symbol = :low, filestring_addl::Vector{String} =
                    Vector{String}(undef, 0), use_highest_posterior_value::Bool = false,
                    input_file_name::String = "")

    if isempty(input_file_name)
        input_file_name = get_forecast_input_file(m, input_type, filestring_addl = filestring_addl)
    end
    println(verbose, :low, "Loading draws from $input_file_name")

    # Load single draw
    if input_type in [:mean, :mode, :mode_draw_shocks]
        if get_setting(m, :sampling_method) == :MH
            params = convert(Vector{Float64}, h5read(input_file_name, "params"))
        elseif get_setting(m, :sampling_method) == :SMC
            if input_type in [:mode, :mode_draw_shocks]
                # If mode is in forecast overrides, want to use the file specified
                if :mode in collect(keys(forecast_input_file_overrides(m)))
                    params = convert(Vector{Float64}, h5read(input_file_name, "params"))
                    # If not, load it from the cloud
                elseif (occursin("smc_paramsmode", input_file_name) &&
                        !use_highest_posterior_value && input_type == :mode)
                    params = convert(Vector{Float64}, h5read(input_file_name, "params"))
                else
                    # Check that input_file_name is correct. Note that if
                    # it already has smc_cloud, then the following code block does not do anything.
                    if occursin("smc_paramsmode", input_file_name) || occursin(".h5", input_file_name)
                        input_file_name = replace(replace(input_file_name,
                                                          "smc_paramsmode" => "smc_cloud"),
                                                  ".h5" => ".jld2")
                        println(verbose, :low, "Switching estimation file of draws to $input_file_name")
                    end

                    cloud = load(input_file_name, "cloud")
                    params = if typeof(cloud) <: Union{DSGE.Cloud,SMC.Cloud}
                        SMC.get_likeliest_particle_value(SMC.Cloud(cloud))
                    else
                        cloud.particles[argmax(get_logpost(cloud))].value
                    end
                end
            else
                error("SMC mean not implemented yet")
            end
        else
            error("Invalid :sampling method specification. Change in setting :sampling_method")
        end
    # Load full distribution
    elseif input_type == :full

        if get_setting(m, :sampling_method) == :MH
            params = map(Float64, h5read(input_file_name, "mhparams"))
        elseif get_setting(m, :sampling_method) == :SMC
            cloud = load(replace(replace(input_file_name, ".h5" => ".jld2"), "smcsave" => "smc_cloud"), "cloud")
            params_unweighted = get_vals(cloud)
            # Re-sample SMC draws according to their weights
            W = load(replace(replace(input_file_name, "smcsave" => "smc_cloud"), "h5" => "jld2"), "W")
            weights = W[:, end]
            inds = SMC.resample(weights)

            params = Matrix{Float64}(params_unweighted[:, inds]')
        else
            error("Invalid :sampling method specification. Change in setting :sampling_method")
        end

    # Load subset of full distribution
    elseif input_type == :subset

        if isempty(subset_inds)
            error("Must supply nonempty range of subset_inds if input_type == :subset")
        else
            if get_setting(m, :sampling_method) == :MH
                params = map(Float64, h5read(input_file_name, "mhparams", (subset_inds, :)))
            elseif get_setting(m, :sampling_method) == :SMC
                params = map(Float64, h5read(input_file_name, "smcparams", (subset_inds, :)))
            else
                throw("Invalid :sampling method specification. Change in setting :sampling_method")
            end
        end
    elseif input_type == :prior
        params = rand(m.parameters, n_forecast_draws(m, :prior))

    elseif input_type == :init || input_type == :init_draw_shocks
        # Return initial parameters of model object

#=        if m.spec == "het_dsge"
            init_parameters!(m, testing_gamma = false)
        else
            init_parameters!(m)
        end =#
        tmp = map(α -> α.value, m.parameters)
        params = convert(Vector{Float64}, tmp)

    end

    return params
end

function load_draws(m::AbstractDSGEModel, input_type::Symbol, block_inds::AbstractRange{Int64};
                    verbose::Symbol = :low,
                    filestring_addl::Vector{String} = Vector{String}(undef, 0),
                    input_file_name::String = "")

    if isempty(input_file_name)
        input_file_name = get_forecast_input_file(m, input_type, filestring_addl = filestring_addl)
    end
    println(verbose, :low, "Loading draws from $input_file_name")

    if input_type in [:full, :subset]
        if isempty(block_inds)
            error("Must supply nonempty range of block_inds for this load_draws method")
        else
            ndraws = length(block_inds)
            params = Vector{Vector{Float64}}(undef, ndraws)
            if get_setting(m, :sampling_method) == :SMC
                params_unweighted = similar(params)
            end
            for (i, j) in zip(1:ndraws, block_inds)
                if get_setting(m, :sampling_method) == :MH
                    params[i] = vec(map(Float64, h5read(input_file_name, "mhparams", (j, :))))
                elseif get_setting(m, :sampling_method) == :SMC

                    params_unweighted[i] = vec(map(Float64, h5read(input_file_name, "smcparams", (j, :))))
                else
                    throw("Invalid :sampling_method setting specification.")
                end
            end

            # Re-sample draws according to weights if the sampling_method was SMC
            if get_setting(m, :sampling_method) == :SMC
                # Re-sample SMC draws according to their weights
                @load replace(replace(input_file_name, "smcsave" => "smc_cloud"), "h5" => "jld2") W
                weights = W[:, end][block_inds]
                inds = SMC.resample(weights)

                params = params_unweighted[inds]
            end

            return params
        end
    elseif input_type == :init_draw_shocks
        ndraws = length(block_inds)
        params = repeat([map(x -> x.value, m.parameters)], ndraws)
    elseif input_type == :mode_draw_shocks
        ndraws = length(block_inds)
        params = repeat([map(x -> x.value, m.parameters)], ndraws)
    elseif input_type == :prior
        ndraws = length(block_inds)
        params = Vector{Vector{Float64}}(undef, ndraws)
        for (i, j) in zip(1:ndraws, block_inds)
            try_again = true
            param_try = vec(rand(m.parameters, ndraws)[:, i])
            # Need to do this to ensure that draws from prior don't return -Inf likelihood
            while try_again
                param_try = vec(rand(m.parameters, ndraws)[:, i])
                try
                    update!(m, param_try)
                    likelihood(m, rand(length(m.observables), 100))
                    try_again = false
                catch
                    try_again = true
                end
            end
            params[i] = param_try
        end
        return params
    else
        error("This load_draws method can only be called with input_type in [:full, :subset]")
    end
end

function load_draws(m::AbstractDSGEVARModel, input_type::Symbol; subset_inds::AbstractRange{Int64} = 1:0,
                    verbose::Symbol = :low,
                    filestring_addl::Vector{String} = Vector{String}(undef, 0),
                    use_highest_posterior_value::Bool = false,
                    input_file_name::String = "")

    if isempty(input_file_name)
        input_file_name = get_forecast_input_file(m, input_type;
                                                  filestring_addl = filestring_addl)
    end

    return load_draws(get_dsge(m), input_type; subset_inds = subset_inds, verbose = verbose,
                      filestring_addl = filestring_addl, use_highest_posterior_value =
                      use_highest_posterior_value,
                      input_file_name = input_file_name)
end

"""
```
forecast_one(m, input_type, cond_type, output_vars; df = DataFrame(),
    subset_inds = 1:0, forecast_string = "", verbose = :low, ...)
```

Compute and save `output_vars` for input draws given by `input_type` and
conditional data case given by `cond_type`.

### Inputs

- `m::AbstractDSGEModel`: model object

- `input_type::Symbol`: one of:

```
  - `:mode`: forecast using the modal parameters only
  - `:mean`: forecast using the mean parameters only
  - `:init`: forecast using the initial parameter values only
  - `:full`: forecast using all parameters (full distribution)
  - `:subset`: forecast using a well-defined user-specified subset of draws
```

- `cond_type::Symbol`: one of:

```
  - `:none`: no conditional data
  - `:semi`: use \"semiconditional data\" - average of quarter-to-date
    observations for high frequency series
  - `:full`: use \"conditional data\" - semiconditional plus nowcasts for
    desired observables
```

- `output_vars::Vector{Symbol}`: vector of desired output variables. See
  `?forecast_one_draw`.

### Keyword Arguments

- `df::DataFrame`: Historical data. If `cond_type in [:semi, :full]`, then the
   final row of `df` should be the period containing conditional data. If not
   provided, will be loaded using `load_data` with the appropriate `cond_type`
- `subset_inds::AbstractRange{Int64}`: indices specifying the draws we want to use. If a
  more sophisticated selection criterion is desired, the user is responsible for
  determining the indices corresponding to that criterion. If `input_type` is
  not `subset`, `subset_inds` will be ignored
- `forecast_string::String`: short string identifying the subset to be
  appended to the output filenames. If `input_type = :subset` and
  `forecast_string` is empty, an error is thrown.
- `only_filter::Bool`: do not run the smoother and only run the filter. This limits the number of
  output variables which can be calculated.
- `verbose::Symbol`: desired frequency of function progress messages printed to
  standard out. One of `:none`, `:low`, or `:high`.
- `check_empty_columns::Bool = true`: check empty columns or not when loading data (if `df` is empty)
- `bdd_fcast::Bool = true`: are we computing the bounded forecasts or not?
- `params::AbstractArray{Float64} = Vector{Float64}(undef, 0)`: parameter draws for the forecast.
     If empty, then we load draws from estimation files implied by the settings in `m`.
- `zlb_method::Symbol`: method for enforcing the zero lower bound. Defaults to `:shock`,
    meaning we use a monetary policy shock to enforce the ZLB.

  Other available methods:
  1. `:temporary_altpolicy` -> use a temporary alternative policy to enforce the ZLB.

- `rerun_smoother::Bool = false`: if true, rerun the conditional forecast when automatically enforcing
    the ZLB as a temporary alternative policy.
- `set_regime_vals_altpolicy::Function`: `Function` that adds new regimes to parameters when
    using temporary alternative policies (if needed). Defaults to identity (which does nothing)
    This function should take as inputs the model object `m` and the total number of regimes
    (after adding the required temporary regimes), i.e. `set_regime_vals_altpolicy(m, n)`. It should then
    set up regime-switching parameters for these new additional regimes.
- `set_info_sets_altpolicy::Function = auto_temp_altpolicy_info_set`: `Function` that automatically updates
    the `tvis_information_set`, e.g. when `zlb_method = :temporary_altpolicy`.
- `show_failed_percent::Bool = false`: prints out the number of failed forecasts, which are returned as NaNs.
    These may occur when the ZLB is not enforced, for example.
- `pegFFR::Bool = false`: peg the nominal FFR at the value specified by `FFRpeg`
- `FFRpeg::Float64 = -0.25/4`: value of the FFR peg
- `H::Int = 4`: number of horizons for which the FFR is pegged
- `testing_carer_kohn::Bool = false`: whether to create a file storing some property of Σ in Carter Kohn

### Outputs

None. Output is saved to files returned by
`get_forecast_output_files(m, input_type, cond_type, output_vars)`.
"""
function forecast_one(m::AbstractDSGEModel{Float64},
                      input_type::Symbol, cond_type::Symbol, output_vars::Vector{Symbol};
                      df::DataFrame = DataFrame(), subset_inds::AbstractRange{Int64} = 1:0,
                      forecast_string::String = "",
                      use_filtered_shocks_in_shockdec::Bool = false,
                      shock_name::Symbol = :none, shock_var_name::Symbol = :none,
                      shock_var_value::Float64 = 0.0, check_empty_columns = true,
                      bdd_fcast::Bool = true, params::AbstractArray{Float64} = Vector{Float64}(undef, 0),
                      zlb_method::Symbol = :shock, set_regime_vals_altpolicy::Function = identity,
                      set_info_sets_altpolicy::Function = auto_temp_altpolicy_info_set,
                      rerun_smoother::Bool = false, catch_smoother_lapack::Bool = false,
                      pegFFR::Bool = false, FFRpeg::Float64 = -0.25/4, H::Int = 4,
                      show_failed_percent::Bool = false, only_filter::Bool = false,
                      verbose::Symbol = :low, testing_carter_kohn::Bool = false,
                      trend_nostates_obs = Array{(0,0)}, trend_nostates_pseudo = Array{(0,0)},
                      full_shock_decomp::Bool = true)

    ### Common Setup

    # Add necessary output_vars and load data
    output_vars, df = prepare_forecast_inputs!(m, input_type, cond_type, output_vars;
                                               df = df, verbose = verbose, bdd_fcast = bdd_fcast,
                                               subset_inds = subset_inds,
                                               check_empty_columns = check_empty_columns, only_filter = only_filter)

    # Get output file names
    forecast_output = Dict{Symbol, Array{Float64}}()
    forecast_output_files = get_forecast_output_files(m, input_type, cond_type, output_vars;
                                                      forecast_string = forecast_string)
    output_dir = rawpath(m, "forecast")

    # Determine if we are regime_switching
    regime_switching = haskey(get_settings(m), :regime_switching) ?
        get_setting(m, :regime_switching) : false
    n_regimes        = regime_switching && haskey(get_settings(m), :n_regimes) ?
        get_setting(m, :n_regimes) : 1

    # Print
    info_print(verbose, :low, "Forecasting input_type = $input_type, cond_type = $cond_type...")
    println(verbose, :low, "Start time: $(now())")
    println(verbose, :low, "Forecast outputs will be saved in $output_dir")


    ### Single-Draw Forecasts

    if input_type in [:mode, :mean, :init]

        elapsed_time = @elapsed let
            if isempty(params)
                params = load_draws(m, input_type; verbose = verbose)
            end
            forecast_output = forecast_one_draw(m, input_type, cond_type, output_vars,
                                                params, df, verbose = verbose,
                                                shock_name = shock_name,
                                                shock_var_name = shock_var_name,
                                                shock_var_value = shock_var_value,
                                                regime_switching = regime_switching,
                                                n_regimes = n_regimes, zlb_method = zlb_method,
                                                set_regime_vals_altpolicy = set_regime_vals_altpolicy,
                                                set_info_sets_altpolicy = set_info_sets_altpolicy,
                                                pegFFR = pegFFR, FFRpeg = FFRpeg, H = H, only_filter = only_filter,
                                                rerun_smoother = rerun_smoother, catch_smoother_lapack = catch_smoother_lapack,
                                                testing_carter_kohn = testing_carter_kohn, trend_nostates_obs = trend_nostates_obs,
                                                trend_nostates_pseudo = trend_nostates_pseudo, full_shock_decomp = full_shock_decomp)

            write_forecast_outputs(m, input_type, output_vars, forecast_output_files,
                                   forecast_output; df = df, block_number = Nullable{Int64}(),
                                   verbose = verbose)
        end

        write_forecast_outputs(m, input_type, output_vars, forecast_output_files,
                               forecast_output; df = df, block_number = Nullable{Int64}(),
                               verbose = verbose)

        total_forecast_time     = elapsed_time
        total_forecast_time_min = total_forecast_time/60
        println(verbose, :low, "\nTotal time to forecast: $total_forecast_time_min minutes")


    ### Multiple-Draw Forecasts

    elseif input_type in [:full, :subset, :prior, :init_draw_shocks, :mode_draw_shocks]
        if get_setting(m, :sampling_method) == :SMC
            subset_cond = input_type == :subset && get_jstep(m, length(subset_inds)) > 1
            other_cond = input_type in [:full, :prior, :init_draw_shocks, :mode_draw_shocks] &&
                get_jstep(m, n_forecast_draws(m, :full)) > 1
            if subset_cond || other_cond
                answer = ""
                invalid_answer = true
                while invalid_answer
                    println("WARNING: The forecast draws are being thinned by $(get_jstep(m, n_forecast_draws(m,:full))). This is not recommended for forecasts using an SMC estimation. Do you want to continue? (y/n)")
                    answer = readline(stdin)
                    if (answer == "y" || answer == "n")
                        invalid_answer = false
                    else
                        println("Invalid answer! Must be `y` or `n`")
                    end
                end

                if answer == "n"
                    throw(error("Forecast halted by user."))
                end
            end
        end

        if !isempty(params) & isa(params, AbstractMatrix)
            m <= Setting(:forecast_ndraws, size(params, 1))
        elseif !isempty(params) & (input_type == :mode_draw_shocks) & (!haskey(get_settings(m), :forecast_ndraws))
            error("If using :mode_draw_shocks with a user-defined matrix of parameter draws, " *
                  "then the user must manually add the Setting :forecast_ndraws " *
                  "to indicate the number of times shocks will be drawn.")
        end

        # Block info
        block_inds, block_inds_thin = forecast_block_inds(m, input_type; subset_inds = subset_inds)
        nblocks = length(block_inds)
        start_block = forecast_start_block(m)

        # Info needed for printing progress
        total_forecast_time = 0.0
        block_verbose = verbose == :none ? :none : :low

        if testing_carter_kohn
            arred = zeros(block_inds[end][end], 246)
        end

        for block = start_block:nblocks
            println(verbose, :low, )
            info_print(verbose, :low, "Forecasting block $block of $nblocks...")
            begin_time = time_ns()

            # Get to work!
            if isempty(params)
                params_for_map = load_draws(m, input_type, block_inds[block]; verbose = verbose)
            elseif input_type == :mode_draw_shocks
                ndraws = length(block_inds[block])
                @assert isa(params, Vector) "To use mode_draw_shocks with params passed in as a keyword, params must be a Vector."
                params_for_map = Vector{Float64}[params for i in block_inds[block]]
            else
                params_for_map = Vector{Float64}[params[i, :] for i in block_inds[block]]
            end

            mapfcn = use_parallel_workers(m) ? pmap : map
            forecast_outputs = mapfcn(param -> forecast_one_draw(m, input_type, cond_type, output_vars,
                                                                 param, df, verbose = verbose,
                                                                 use_filtered_shocks_in_shockdec =
                                                                 use_filtered_shocks_in_shockdec,
                                                                 shock_name = shock_name,
                                                                 shock_var_name = shock_var_name,
                                                                 shock_var_value = shock_var_value,
                                                                 regime_switching = regime_switching,
                                                                 n_regimes = n_regimes, zlb_method = zlb_method,
                                                                 set_regime_vals_altpolicy = set_regime_vals_altpolicy,
                                                                 set_info_sets_altpolicy = set_info_sets_altpolicy,
                                                                 pegFFR = pegFFR, FFRpeg = FFRpeg, H = H, only_filter = only_filter,
                                                                 rerun_smoother = rerun_smoother,
                                                                 catch_smoother_lapack = catch_smoother_lapack,
                                                                 testing_carter_kohn = testing_carter_kohn),
                                      params_for_map)

            # Assemble outputs from this block and write to file
            forecast_outputs = convert(Vector{Dict{Symbol, Array{Float64}}}, forecast_outputs)

            if testing_carter_kohn && get_setting(m, :forecast_smoother) == :carter_kohn
                for i in 1:length(forecast_outputs)
                    arred[block_inds[block][i],1:length(forecast_outputs[i][:conded])] = forecast_outputs[i][:conded]
                    delete!(forecast_outputs[i], :conded)
                end
            end


            forecast_output = assemble_block_outputs(forecast_outputs; show_failed_percent = show_failed_percent)
            #@show forecast_output.keys
            write_forecast_outputs(m, input_type, output_vars, forecast_output_files,
                                   forecast_output; df = df, block_number = Nullable(block),
                                   verbose = block_verbose, block_inds = block_inds_thin[block],
                                   subset_inds = subset_inds)
            GC.gc()

            # Calculate time to complete this block, average block time, and
            # expected time to completion
            block_time = (time_ns() - begin_time)/1e9
            total_forecast_time += block_time
            total_forecast_time_min     = total_forecast_time/60
            blocks_elapsed              = block - start_block + 1
            expected_time_remaining     = (total_forecast_time/blocks_elapsed)*(nblocks - block)
            expected_time_remaining_min = expected_time_remaining/60

            println(verbose, :low, "\nCompleted $block of $nblocks blocks.")
            println(verbose, :low, "Total time elapsed: $total_forecast_time_min minutes")
            println(verbose, :low, "Expected time remaining: $expected_time_remaining_min minutes")
        end # of loop through blocks

        if testing_carter_kohn && get_setting(m, :forecast_smoother) == :carter_kohn
            jldopen("conded_svd_diagm.jld2", true, true, true, IOStream) do file
                write(file, "conds", arred)
            end
        end
    end # of input_type
    combine_raw_forecast_output_and_metadata(m, forecast_output_files, verbose = verbose)

    println(verbose, :low, "\nForecast complete: $(now())")
end

"""
```
forecast_one_draw(m, input_type, cond_type, output_vars, params, df;
    verbose = :low, only_filter = false)
```

Compute `output_vars` for a single parameter draw, `params`. Called by
`forecast_one`.

### Inputs

- `m::AbstractDSGEModel{Float64}`: model object
- `input_type::Symbol`: See `?forecast_one`.
- `cond_type::Symbol`: See `?forecast_one`.
- `output_vars::Vector{Symbol}`: vector of desired output variables. See Outputs
  section
- `params::Vector{Float64}`: parameter vector
- `df::DataFrame`: historical data.
- `verbose::Symbol`: desired frequency of function progress messages printed to
  standard out. One of `:none`, `:low`, or `:high`.

### Output

- `forecast_outputs::Dict{Symbol, Array{Float64}}`: dictionary of forecast
  outputs. Keys are `output_vars`, which is some subset of:

```
  - `:histstates`: `Matrix{Float64}` of smoothed historical states
  - `:histobs`: `Matrix{Float64}` of smoothed historical data
  - `:histpseudo`: `Matrix{Float64}` of smoothed historical
    pseudo-observables
  - `:histshocks`: `Matrix{Float64}` of smoothed historical shocks
  - `:forecaststates`: `Matrix{Float64}` of forecasted states
  - `:forecastobs`: `Matrix{Float64}` of forecasted observables
  - `:forecastpseudo`: `Matrix{Float64}` of forecasted pseudo-observables
  - `:forecastshocks`: `Matrix{Float64}` of forecasted shocks
  - `:bddforecaststates`, `:bddforecastobs`, `:bddforecastpseudo`, and
    `:bddforecastshocks`: `Matrix{Float64}`s of forecasts where we enforce
    the zero lower bound to be `forecast_zlb_value(m)`
  - `:shockdecstates`: `Array{Float64, 3}` of state shock decompositions
  - `:shockdecobs`: `Array{Float64, 3}` of observable shock decompositions
  - `:shockdecpseudo`: `Array{Float64, 3}` of pseudo-observable shock
    decompositions
  - `:dettrendstates`: `Matrix{Float64}` of state deterministic trends
  - `:dettrendobs`: `Matrix{Float64}` of observable deterministic trends
  - `:dettrendpseudo`: `Matrix{Float64}` of pseudo-observable deterministic
    trends
  - `:trendstates`: `Vector{Float64}` of state trends, i.e. the `CCC` vector
  - `:trendobs`: `Vector{Float64}` of observable trends, i.e. the `DD` vector
  - `:trendpseudo`: `Vector{Float64}` of pseudo-observable trends, i.e. the
    `DD_pseudo` vector
  - `:irfstates`: `Array{Float64, 3}` of state impulse responses
  - `:irfobs`: `Array{Float64, 3}` of observable impulse responses
  - `:irfpseudo`: `Array{Float64, 3}` of pseudo-observable impulse responses
```
"""
function forecast_one_draw(m::AbstractDSGEModel{Float64}, input_type::Symbol, cond_type::Symbol,
                           output_vars::Vector{Symbol}, params::Vector{Float64}, df::DataFrame; verbose::Symbol = :low,
                           use_filtered_shocks_in_shockdec::Bool = false,
                           shock_name::Symbol = :none, shock_var_name::Symbol = :none,
                           shock_var_value::Float64 = 0.0, zlb_method::Symbol = :shock,
                           set_regime_vals_altpolicy::Function = identity,
                           set_info_sets_altpolicy::Function = auto_temp_altpolicy_info_set,
                           pegFFR::Bool = false, FFRpeg::Float64 = -0.25/4, H::Int = 4,
                           regime_switching::Bool = false, n_regimes::Int = 1, only_filter::Bool = false,
                           filter_smooth::Bool = false, rerun_smoother::Bool = false,
                           catch_smoother_lapack::Bool = false, testing_carter_kohn::Bool = false,
                           return_loglh::Bool = false, trend_nostates_obs = Array{(0,0)},
                           trend_nostates_pseudo = Array{(0,0)}, full_shock_decomp::Bool = false)
    ### Setup

    # Re-initialize model indices if forecasting under an alternative policy
    # rule (in case new states or equations were added)
    if alternative_policy(m).key != :historical
        init_model_indices!(m)
    end

    # Time-Varying Information Set?
    tvis = haskey(get_settings(m), :tvis_information_set) ? !isempty(get_setting(m, :tvis_information_set)) : false

    # Are we only running IRFs?
    output_prods = map(get_product, output_vars)
    irfs_only = all(x -> x == :irf, output_prods)

    # Compute state space
    update!(m, params) # Note that params is a Vector{Float64}, not a ParameterVector. This `update!` infers if the forecast is regime-switching if length(params) > length(m.parameters)

    system = compute_system(m; tvis = tvis)

    # Initialize output dictionary
    forecast_output = Dict{Symbol, Array{Float64}}()

    # Decide whether to draw states/shocks in smoother/forecast
    uncertainty_override = forecast_uncertainty_override(m)
    uncertainty = if Nullables.isnull(uncertainty_override)
        if input_type in [:init, :mode, :mean, :prior]
            false
        elseif input_type in [:full, :subset, :mode_draw_shocks, :init_draw_shocks]
            true
        end
    else
        get(uncertainty_override)
    end

    ### 1. Smoothed Histories

    # Must run smoother for conditional data in addition to explicit cases
    hist_vars = [:histstates, :histpseudo, :histshocks, :histstdshocks]
    shockdec_vars = [:shockdecstates, :shockdecpseudo, :shockdecobs]
    dettrend_vars = [:dettrendstates, :dettrendpseudo, :dettrendobs]
    smooth_vars = vcat(hist_vars, shockdec_vars, dettrend_vars)
    hists_to_compute = intersect(output_vars, smooth_vars)

    run_smoother = (!isempty(hists_to_compute) ||
        (cond_type in [:semi, :full] && !irfs_only)) && !only_filter
    if run_smoother
        # Call smoother
        if filter_smooth && (get_setting(m, :forecast_smoother) == :carter_kohn)
            kal = filter(m, df, system; cond_type = cond_type)
            histstates, histshocks, histpseudo, initial_states =
                smooth(m, df, system; cond_type = cond_type, draw_states = uncertainty,
                       s_pred = kal[:s_pred], P_pred = kal[:P_pred], s_filt = kal[:s_filt], P_filt = kal[:P_filt])

        else
            histstates, histshocks, histpseudo, initial_states =
                smooth(m, df, system; cond_type = cond_type, draw_states = uncertainty)
        end

        # For conditional data, transplant the obs/state/pseudo vectors from hist to forecast
        if cond_type in [:full, :semi]
            T = n_mainsample_periods(m)

            forecast_output[:histstates] = transplant_history(histstates, T)
            forecast_output[:histshocks] = transplant_history(histshocks, T)
            forecast_output[:histpseudo] = transplant_history(histpseudo, T)
        else
            forecast_output[:histstates] = histstates
            forecast_output[:histshocks] = histshocks
            forecast_output[:histpseudo] = histpseudo
        end

        # Standardize shocks if desired
        if :histstdshocks in output_vars
            if regime_switching
                start_date = max(date_mainsample_start(m), df[1, :date]) # use mainsample b/c shocks only includes mainsample, no presample
                end_date   = prev_quarter(date_forecast_start(m))
                regime_inds = regime_indices(m, start_date, end_date)
                if regime_inds[1][1] < 1
                    regime_inds[1] = 1:regime_inds[1][end]
                end
                forecast_output[:histstdshocks] =
                    standardize_shocks(forecast_output[:histshocks],
                                       Matrix{eltype(system[1, :QQ])}[system[i, :QQ] for i in 1:n_regimes], regime_inds)
            else
                forecast_output[:histstdshocks] = standardize_shocks(forecast_output[:histshocks], system[:QQ])
            end
        end
    end

    ### 2. Forecasts

    unbddforecast_vars = [:forecaststates, :forecastobs, :forecastpseudo, :forecastshocks, :forecaststdshocks]
    bddforecast_vars = [:bddforecaststates, :bddforecastobs, :bddforecastpseudo, :bddforecastshocks, :bddforecaststdshocks]
    forecast_vars = vcat(unbddforecast_vars, bddforecast_vars)
    forecasts_to_compute = intersect(output_vars, forecast_vars)

    if !isempty(forecasts_to_compute)
        # Get initial forecast state vector s_T
        s_T = if run_smoother && !return_loglh # ONLY THIS BRANCH WORKS FOR REGIME SWITCHING
            # The last smoothed state is either s_{T|T} (if !uncertainty) or
            # drawn from N(s_{T|T}, P_{T|T}) (if uncertainty)
            histstates[:, end]
        else
            kal = Kalman(Vector{Float64}(undef,0), Matrix{Float64}(undef, 0, 0), Array{Float64}(undef, 0, 0, 0), Matrix{Float64}(undef, 0, 0), Array{Float64}(undef, 0, 0, 0), Vector{Float64}(undef, 0), Array{Float64}(undef, 0, 0, 0), Vector{Float64}(undef, 0), Array{Float64}(undef, 0, 0, 0))
            try
                kal = filter(m, df, system; cond_type = cond_type)
            catch err
                if isa(err, DomainError)
                    return Dict{Symbol, Array{Float64}}()
                else
                    rethrow(err)
                end
            end
            if uncertainty
                # If we want to draw s_T but haven't run the smoother, draw from
                # N(s_{T|T}, P_{T|T}) directly
                U, singular_values, _ = svd(kal[:P_T])
                dist = DegenerateMvNormal(kal[:s_T], U*Matrix(Diagonal(sqrt.(singular_values))))
                rand(dist)
            else
                # If we don't want to draw s_T, simply use the mean s_{T|T}
                kal[:s_T]
            end
        end

        # Re-solve model with alternative policy rule, if applicable
        apply_altpolicy = alternative_policy(m).solve != solve
        if apply_altpolicy

            if haskey(get_settings(m), :skip_altpolicy_state_init) ? !get_setting(m, :skip_altpolicy_state_init) : true
                # Adjust the initial state vector for pgap and ygap
                if haskey(m.settings, :pgap_type) && haskey(get_settings(m), :pgap_value)
                    update_hist = true
                    _pgap_type = get_setting(m, :pgap_type)
                    if _pgap_type == :ngdp
                        _, s_T = ngdp_forecast_init(m, zeros(0, 0), s_T, cond_type = cond_type)
                    elseif _pgap_type == :ait
                        _, s_T = ait_forecast_init(m, zeros(0, 0), s_T, cond_type = cond_type)
                    else
                        update_hist = false
                    end
                    if update_hist # Don't do this computation unless pgap_type in [:ngdp, :ait]
                        histstates[:, end] = s_T
                        ZZ_pseudo = isa(system, RegimeSwitchingSystem) ? system[1, :ZZ_pseudo] : system[:ZZ_pseudo] # currently assumes these are
                        DD_pseudo = isa(system, RegimeSwitchingSystem) ? system[1, :DD_pseudo] : system[:DD_pseudo] # the same across regimes
                        histpseudo[:, end] = ZZ_pseudo * s_T + DD_pseudo
                    end
                end

                if haskey(m.settings, :ygap_type) && haskey(get_settings(m), :ygap_value) &&
                    haskey(m.settings, :pgap_type) && haskey(get_settings(m), :pgap_value)
                    update_hist = true
                    _ygap_type = get_setting(m, :ygap_type)
                    _pgap_type = get_setting(m, :pgap_type)
                    if _ygap_type == _pgap_type == :smooth_ait_gdp
                        _, s_T = smooth_ait_gdp_forecast_init(m, zeros(0, 0), s_T, cond_type = cond_type)
                    elseif _ygap_type == _pgap_type == :smooth_ait_gdp_alt
                        _, s_T = smooth_ait_gdp_alt_forecast_init(m, zeros(0, 0), s_T, cond_type = cond_type)
                    elseif (_ygap_type == _pgap_type == :flexible_ait) &&
                        (haskey(m.settings, :flexible_ait_policy_change) ?
                         !get_setting(m, :flexible_ait_policy_change) : true) &&
                         (haskey(m.settings, :initialize_pgap_ygap) ? get_setting(m, :initialize_pgap_ygap) : false)
                        _, s_T = flexible_ait_forecast_init(m, zeros(0, 0), s_T, cond_type = cond_type)
                    else
                        update_hist = false
                    end
                    if update_hist # Don't do this computation unless pgap_type and ygap_type are okay
                        histstates[:, end] = s_T
                        ZZ_pseudo = isa(system, RegimeSwitchingSystem) ? system[1, :ZZ_pseudo] : system[:ZZ_pseudo] # currently assumes these are
                        DD_pseudo = isa(system, RegimeSwitchingSystem) ? system[1, :DD_pseudo] : system[:DD_pseudo] # the same across regimes
                        histpseudo[:, end] = ZZ_pseudo * s_T + DD_pseudo
                    end
                end

                if haskey(m.settings, :Rref_type)
                    if get_setting(m, :ygap_type) == :rw && get_setting(m, :pgap_type) == :rw &&
                        get_setting(m, :Rref_type) == :rw
                        _, s_T = rw_forecast_init(m, zeros(0, 0), s_T, cond_type = cond_type)
                        histstates[:, end] = s_T
                        ZZ_pseudo = isa(system, RegimeSwitchingSystem) ? system[1, :ZZ_pseudo] : system[:ZZ_pseudo] # currently assumes these are
                        DD_pseudo = isa(system, RegimeSwitchingSystem) ? system[1, :DD_pseudo] : system[:DD_pseudo] # the same across regimes
                        histpseudo[:, end] = ZZ_pseudo * s_T + DD_pseudo
                    end
                end
            end
        end

        # 2A. Unbounded forecasts
        if !isempty(intersect(output_vars, unbddforecast_vars))
            if pegFFR
                nshocks = size(system[:RRR], 2)
                nstates = size(system[:TTT], 1)
                shocks = m.exogenous_shocks
                PsiR1 = 0
                PsiR2 = zeros(nstates)
                PsiR2[m.endogenous_states[:R_t]] = 1.
                Rht = system[:RRR][:, vcat(shocks[:rm_sh], shocks[:rm_shl1]:shocks[Symbol("rm_shl$H")])]
                bb = zeros(H+1, 1)
                MH = zeros(H+1, H+1)
                for hh = 1:H+1
                    bb[hh, 1] = (FFRpeg - PsiR1 - PsiR2'*(system[:TTT])^hh*s_T)[1]
                    MH[hh, :] = PsiR2'*(system[:TTT])^(hh-1)*Rht
                end
                monshocks = MH\bb
                etpeg = zeros(nshocks, forecast_horizons(m))
                etpeg[vcat(shocks[:rm_sh], shocks[:rm_shl1]:shocks[Symbol("rm_shl$H")]), 1] = monshocks
                forecaststates, forecastobs, forecastpseudo, forecastshocks =
                    forecast(system, s_T, etpeg; cond_type = cond_type, enforce_zlb = false, draw_shocks = uncertainty)
                println("The forecasted interest rate path is $(forecastobs[m.observables[:obs_nominalrate], :])")
            else
                forecaststates, forecastobs, forecastpseudo, forecastshocks =
                    forecast(m, system, s_T;
                             cond_type = cond_type, enforce_zlb = false, draw_shocks = uncertainty)
            end

            # For conditional data, transplant the obs/state/pseudo vectors from hist to forecast
            # NOTE: ZZ REGIME SWITCHING NOT FULLY SUPPORTED, SO JUST TAKE THE LAST ZZ IN THE SYSTEM
            #       OR THE ZZ SPECIFIED BY THE FOLLOWING SETTING
            cond_ZZ_regime = if regime_switching
                if haskey(get_settings(m), :cond_ZZ_regime)
                    get_setting(m, :cond_ZZ_regime)
                else # an educated guess at the correct conditional regime
                    if only_filter
                        @warn "The setting :cond_ZZ_regime is not specified. Making an educated guess..."
                    end
                    if get_setting(m, :n_cond_regimes) > 0 # if regime-switching during conditional horizon, need to add them and
                        get_setting(m, :reg_forecast_start) + get_setting(m, :n_cond_regimes) - 1 # minus 1 for reg_forecast_start
                    else # otherwise no regime-switching during conditional horizon, so just use the same regime as forecast start date
                        get_setting(m, :reg_forecast_start)
                    end
                end
            else
                1
            end
            if cond_type in [:full, :semi] && !only_filter
                forecast_output[:forecaststates] = transplant_forecast(histstates, forecaststates, T)
                forecast_output[:forecastshocks] = transplant_forecast(histshocks, forecastshocks, T)
                forecast_output[:forecastpseudo] = transplant_forecast(histpseudo, forecastpseudo, T)
                forecast_output[:forecastobs]    =
                    transplant_forecast_observables(histstates, forecastobs,
                                                    regime_switching ? system[cond_ZZ_regime] : system, T)
            elseif cond_type in [:full, :semi] && only_filter
                if regime_switching
                    cond_obs                         = system[cond_ZZ_regime, :ZZ] * s_T + system[cond_ZZ_regime, :DD]
                    cond_pseudo                      = system[cond_ZZ_regime, :ZZ_pseudo] * s_T + system[cond_ZZ_regime, :DD_pseudo]
                    forecast_output[:forecaststates] = hcat(s_T,         forecaststates)
                    forecast_output[:forecastobs]    = hcat(cond_obs,    forecastobs)
                    forecast_output[:forecastpseudo] = hcat(cond_pseudo, forecastpseudo)
                else
                    cond_obs                         = system[:ZZ] * s_T + system[:DD]
                    cond_pseudo                      = system[:ZZ_pseudo] * s_T + system[:DD_pseudo]
                    forecast_output[:forecaststates] = hcat(s_T,         forecaststates)
                    forecast_output[:forecastobs]    = hcat(cond_obs,    forecastobs)
                    forecast_output[:forecastpseudo] = hcat(cond_pseudo, forecastpseudo)
                    # forecastshocks cannot be calculated w/out smoother
                end
            else
                forecast_output[:forecaststates] = forecaststates
                forecast_output[:forecastshocks] = forecastshocks
                forecast_output[:forecastpseudo] = forecastpseudo
                forecast_output[:forecastobs]    = forecastobs
            end

            # Standardize shocks if desired
            if :forecaststdshocks in output_vars
                forecast_output[:forecaststdshocks] =
                    standardize_shocks(forecast_output[:forecastshocks],
                                       regime_switching ? system[cond_ZZ_regime, :QQ] : system[:QQ])
            end
        end

        # 2B. Bounded forecasts

        if !isempty(intersect(output_vars, bddforecast_vars))
            if pegFFR
                nshocks = size(system[:RRR], 2)
                nstates = size(system[:TTT], 1)
                shocks = m.exogenous_shocks
                PsiR1 = 0
                PsiR2 = zeros(nstates)
                PsiR2[m.endogenous_states[:R_t]] = 1.
                Rht = system[:RRR][:, vcat(shocks[:rm_sh], shocks[:rm_shl1]:shocks[Symbol("rm_shl$H")])]
                bb = zeros(H+1, 1)
                MH = zeros(H+1, H+1)
                @show s_T[m.endogenous_states[:R_t]] + 400*log(m[:Rstarn])
                for hh = 1:H+1
                    bb[hh, 1] = (FFRpeg - PsiR1 - PsiR2'*(system[:TTT])^hh*s_T)[1]
                    MH[hh, :] = PsiR2'*(system[:TTT])^(hh-1)*Rht
                end
                monshocks = MH\bb
                etpeg = zeros(nshocks, forecast_horizons(m))
                etpeg[vcat(shocks[:rm_sh], shocks[:rm_shl1]:shocks[Symbol("rm_shl$H")]), 1] = monshocks
                forecaststates, forecastobs, forecastpseudo, forecastshocks = forecast(system, s_T, etpeg)
                @show forecaststates[m.endogenous_states[:R_t], :]
                @show forecastobs[m.observables[:obs_nominalrate], :]
            elseif zlb_method == :temporary_altpolicy
                altpolicy = alternative_policy(m).key
                @assert altpolicy in [:historical, :ait, :ngdp, :smooth_ait_gdp, :smooth_ait_gdp_alt,
                                      :flexible_ait] "altpolicy must be permanent " *
                    "and among [:historical, :ait, :ngdp, :smooth_ait_gdp, :smooth_ait_gdp_alt] for this method of enforcing the ZLB."

                # Run the unbounded forecast if they haven't already been computed
                if isempty(intersect(output_vars, unbddforecast_vars))
                    forecaststates, forecastobs, forecastpseudo, forecastshocks =
                        forecast(m, system, s_T;
                                 cond_type = cond_type, enforce_zlb = false, draw_shocks = uncertainty)
                end

                # Now run the ZLB enforcing forecast
                if rerun_smoother
                    zlb_enforced_output =
                        forecast(m, altpolicy, s_T, forecaststates, forecastobs, forecastpseudo, forecastshocks;
                                 cond_type = cond_type,
                                 set_zlb_regime_vals = set_regime_vals_altpolicy,
                                 set_info_sets_altpolicy = set_info_sets_altpolicy,
                                 rerun_smoother = rerun_smoother, df = df,
                                 draw_states = uncertainty,
                                 histstates = histstates, histshocks = histshocks, histpseudo = histpseudo,
                                 initial_states = initial_states)
                    forecaststates, forecastobs, forecastpseudo, histstates, histshocks, histpseudo, initial_states =
                        zlb_enforced_output

                    if cond_type in [:full, :semi]
                        T = n_mainsample_periods(m)

                        forecast_output[:histstates] = transplant_history(histstates, T)
                        forecast_output[:histshocks] = transplant_history(histshocks, T)
                        forecast_output[:histpseudo] = transplant_history(histpseudo, T)
                    else
                        forecast_output[:histstates] = histstates
                        forecast_output[:histshocks] = histshocks
                        forecast_output[:histpseudo] = histpseudo
                    end

                    # Standardize shocks if desired
                    if :histstdshocks in output_vars
                        if regime_switching
                            start_date = max(date_mainsample_start(m), df[1, :date]) # use mainsample b/c shocks only includes mainsample, no presample
                            end_date   = prev_quarter(date_forecast_start(m))
                            regime_inds = regime_indices(m, start_date, end_date)
                            if regime_inds[1][1] < 1
                                regime_inds[1] = 1:regime_inds[1][end]
                            end
                            forecast_output[:histstdshocks] =
                                standardize_shocks(forecast_output[:histshocks],
                                                   Matrix{eltype(system[1, :QQ])}[system[i, :QQ] for i in 1:n_regimes], regime_inds)
                        else
                            forecast_output[:histstdshocks] = standardize_shocks(forecast_output[:histshocks], system[:QQ])
                        end
                    end
                else
                    forecaststates, forecastobs, forecastpseudo =
                        forecast(m, altpolicy, s_T, forecaststates, forecastobs, forecastpseudo, forecastshocks;
                                 cond_type = cond_type,
                                 set_zlb_regime_vals = set_regime_vals_altpolicy,
                                 set_info_sets_altpolicy = set_info_sets_altpolicy)
                end
            else
                forecaststates, forecastobs, forecastpseudo, forecastshocks =
                    forecast(m, system, s_T;
                             cond_type = cond_type, enforce_zlb = true, draw_shocks = uncertainty)
            end

            # For conditional data, transplant the obs/state/pseudo vectors from hist to forecast
            # NOTE: ZZ REGIME SWITCHING NOT FULLY SUPPORTED, SO JUST TAKE THE LAST ZZ IN THE SYSTEM
            #       OR THE ZZ SPECIFIED BY THE FOLLOWING SETTING
            cond_ZZ_regime = if regime_switching
                if haskey(get_settings(m), :cond_ZZ_regime)
                    get_setting(m, :cond_ZZ_regime)
                else # an educated guess at the correct conditional regime
                    if only_filter
                        @warn "The setting :cond_ZZ_regime is not specified. Making an educated guess..."
                    end
                    if get_setting(m, :n_cond_regimes) > 0 # if regime-switching during conditional horizon, need to add them and
                        get_setting(m, :reg_forecast_start) + get_setting(m, :n_cond_regimes) - 1 # minus 1 for reg_forecast_start
                    else # otherwise no regime-switching during conditional horizon, so just use the same regime as forecast start date
                        get_setting(m, :reg_forecast_start)
                    end
                end
            else
                1
            end
            if cond_type in [:full, :semi] && !only_filter
                forecast_output[:bddforecaststates] = transplant_forecast(histstates, forecaststates, T)
                forecast_output[:bddforecastshocks] = transplant_forecast(histshocks, forecastshocks, T)
                forecast_output[:bddforecastpseudo] = transplant_forecast(histpseudo, forecastpseudo, T)
                forecast_output[:bddforecastobs]    =
                    transplant_forecast_observables(histstates, forecastobs,
                                                    regime_switching ? system[cond_ZZ_regime] : system, T)
            elseif cond_type in [:full, :semi] && only_filter
                forecast_output[:bddforecaststates] = hcat(s_T,forecaststates)
                #forecast_output[:bddforecastshocks] = hcat(s_T,forecastshocks)
                #forecast_output[:bddforecastpseudo] = hcat(s_T,forecastpseudo)
                cond_var = system[cond_ZZ_regime][:ZZ]*s_T .+ system[cond_ZZ_regime][:DD]
                forecast_output[:bddforecastobs]    = hcat(cond_var,forecastobs)
            else
                forecast_output[:bddforecaststates] = forecaststates
                forecast_output[:bddforecastshocks] = forecastshocks
                forecast_output[:bddforecastpseudo] = forecastpseudo
                forecast_output[:bddforecastobs]    = forecastobs
            end

            # Standardize shocks if desired
            if :bddforecaststdshocks in output_vars
                forecast_output[:bddforecaststdshocks] =
                    standardize_shocks(forecast_output[:bddforecastshocks],
                                       regime_switching ? system[cond_ZZ_regime, :QQ] : system[:QQ])
            end
        end
    end

    if return_loglh
        return kal[:loglh]
    end

    ### 3. Shock Decompositions

    shockdecs_to_compute = intersect(output_vars, shockdec_vars)

    if !isempty(shockdecs_to_compute)
        histshocks_shockdec = if use_filtered_shocks_in_shockdec
            filter_shocks(m, df, system, cond_type = cond_type)
        else
            histshocks
        end

        start_date = max(date_mainsample_start(m), df[1, :date]) # smooth doesn't return presample
        end_date   = if cond_type in [:semi, :full] # end date of histshocks includes conditional periods
            max(date_conditional_end(m), prev_quarter(date_forecast_start(m)))
        else
            prev_quarter(date_forecast_start(m)) # this is the end date of history period
        end

        if haskey(m.settings, :old_shock_decs) && get_setting(m, :old_shock_decs)
            # Must be using TV cred system
            old_system = 0
            for i in sort!(collect(keys(get_setting(m, :regime_eqcond_info))))
                if get_setting(m, :regime_eqcond_info)[i].alternative_policy.key == :zero_rate
                    old_system = regime_indices(m, start_date, end_date)[i][1]
                    break
                end
            end
            shockdecstates, shockdecobs, shockdecpseudo = shock_decompositions(m, system, old_system, histshocks_shockdec, start_date, end_date, cond_type, full_shock_decomp = full_shock_decomp)
        else
            shockdecstates, shockdecobs, shockdecpseudo = isa(system, RegimeSwitchingSystem) ?
                shock_decompositions(m, system, histshocks_shockdec, start_date, end_date, cond_type) :
                shock_decompositions(m, system, histshocks_shockdec)
        end

        forecast_output[:shockdecstates] = shockdecstates
        forecast_output[:shockdecobs]    = shockdecobs
        forecast_output[:shockdecpseudo] = shockdecpseudo
    end


    ### 4. Trend

    trend_vars = [:trendstates, :trendobs, :trendpseudo]
    trends_to_compute = intersect(output_vars, trend_vars)

    if !isempty(trends_to_compute)
        trendstates, trendobs, trendpseudo = if regime_switching &&
            (haskey(get_settings(m), :time_varying_trends) ? get_setting(m, :time_varying_trends) : false)

            start_date = max(date_mainsample_start(m), df[1, :date]) # smooth doesn't return presample
            end_date   = if cond_type in [:semi, :full] # end date of histshocks includes conditional periods
                max(date_conditional_end(m), prev_quarter(date_forecast_start(m)), date_forecast_end(m))
            else
                max(prev_quarter(date_forecast_start(m)), date_forecast_end(m)) # this is the end date of history period
            end # need to add date_forecast_end to end_date to avoid dimension mismatch in write_forecast_output

            trends(m, system, start_date, end_date, cond_type)
        else
            trends(system)
        end

        forecast_output[:trendstates] = trendstates
        forecast_output[:trendobs]    = trendobs
        forecast_output[:trendpseudo] = trendpseudo

        if trend_nostates_obs != Array{(0,0)}
            forecast_output[:trendobs] = trend_nostates_obs
        end
        if trend_nostates_pseudo != Array{(0,0)}
            forecast_output[:trendpseudo] = trend_nostates_pseudo
        end
    end


    ### 5. Deterministic Trend

    dettrends_to_compute = intersect(output_vars, dettrend_vars)

    if !isempty(dettrends_to_compute)
        if isa(system, RegimeSwitchingSystem)
            start_date = max(date_mainsample_start(m), df[1, :date]) # smooth doesn't return presample
            end_date   = if cond_type in [:semi, :full] # end date of histshocks includes conditional periods
                max(date_conditional_end(m), prev_quarter(date_forecast_start(m)))
            else
                prev_quarter(date_forecast_start(m)) # this is the end date of history period
            end

            if haskey(m.settings, :old_shock_decs) && get_setting(m, :old_shock_decs) && full_shock_decomp
                # Must be using TV cred system
                old_system = 0
                for i in sort!(collect(keys(get_setting(m, :regime_eqcond_info))))
                    # TODO: Generalize beyond zero_rate to any temporary or permant alternative policy
                    if get_setting(m, :regime_eqcond_info)[i].alternative_policy.key == :zero_rate
                        old_system = regime_indices(m, start_date, end_date)[i][1]
                        break
                    end
                end
                dettrendstates, dettrendobs, dettrendpseudo =
                    deterministic_trends(m, system, old_system, initial_states, start_date, end_date, cond_type)
            else
                dettrendstates, dettrendobs, dettrendpseudo =
                    deterministic_trends(m, system, initial_states, start_date, end_date, cond_type)
            end
        else
            dettrendstates, dettrendobs, dettrendpseudo = deterministic_trends(m, system, initial_states)
        end

        forecast_output[:dettrendstates] = dettrendstates
        forecast_output[:dettrendobs]    = dettrendobs
        forecast_output[:dettrendpseudo] = dettrendpseudo
    end


    ### 6. Impulse Responses

    irf_vars = [:irfstates, :irfobs, :irfpseudo]
    irfs_to_compute = intersect(output_vars, irf_vars)

    if !isempty(irfs_to_compute)
        if shock_name!=:none
            irfstates, irfobs, irfpseudo = impulse_responses(m, system, impulse_response_horizons(m),
                                                             shock_name,
                                                             shock_var_name,
                                                             shock_var_value)
        else
            irfstates, irfobs, irfpseudo = impulse_responses(m, system)
        end
        forecast_output[:irfstates] = irfstates
        forecast_output[:irfobs] = irfobs
        forecast_output[:irfpseudo] = irfpseudo
    end

    ### Return only desired output_vars

    for key in keys(forecast_output)
        if !(key in output_vars)
            delete!(forecast_output, key)
        end
    end

    if testing_carter_kohn && input_type == :full && get_setting(m, :forecast_smoother) == :carter_kohn
        forecast_output[:conded] = conded
    end
    return forecast_output
end
