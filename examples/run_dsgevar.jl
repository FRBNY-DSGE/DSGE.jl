using DSGE, ModelConstructors, Random, Plots, Distributed
using Plots.PlotMeasures

# What do you want to do?
estimate_dsgevar  = true   # Estimate a DSGEVAR using SMC
use_estim_output  = true   # Use posterior from SMC estimation for following code. Otherwise, use default parameters.
get_VAR_system    = true   # Compute VAR coefficients and innovations-covariance matrix
do_modal_irf      = true   # Compute IRFs using modal parameters
compare_modal_irf = true   # Plot the modal DSGEVAR λ = ∞ rotation IRF and the actual DSGE IRF for observables
do_full_band_irf  = true   # Compute IRFs using parameters drawn from a distribution
create_meansbands = false  # Save full_band_irfs to MeansBands
do_parallel       = false  # Use parallel workers
n_workers         = 10
dsgevar_λ         = 0.5    # What λ do you want to use?

Random.seed!(20201793)

# Set up DSGE and load data
m = AnSchorfheide("ss0")
m[:e_R] = 0. # DSGEVAR does not use
m[:e_π] = 0. # measurement error with
m[:e_y] = 0. # subspec ss3
modal_paras = map(x -> x.value, m.parameters)
forecast_string = "test_dsgevar" # Change this to an empty string if you don't want an identifier for saved output
m <= Setting(:sampling_method, :SMC)
m <= Setting(:impulse_response_horizons, 20)

m <= Setting(:n_particles, 50)          # a small number just to make the estimation and IRFs faster.
m <= Setting(:use_fixed_schedule, true) # Fix schedule to fewer iterations
m <= Setting(:n_Φ, 50)                  # to speed up the estimation further

m <= Setting(:date_mainsample_start, DSGE.quartertodate("1960-Q1")) # Set bounds on the periods of
m <= Setting(:date_presample_start, DSGE.quartertodate("1959-Q3"))  # data being loaded so we don't
m <= Setting(:date_forecast_start, DSGE.quartertodate("2019-Q3"))   # get NaNs.
data = df_to_matrix(m, load_data(m))    # the rows of the data MUST correspond
                                        # to the exact order of the observables specified by the DSGEVAR.

# Construct a DSGEVAR from AnSchorfheide
dsgevar = DSGEVAR(m, collect(keys(m.exogenous_shocks)), "ss3")
DSGE.update!(dsgevar, λ = dsgevar_λ)

#= ss3 prespecifies the observables, lags, and λ
Equivalently, we could have run

```
dsgevar = DSGEVAR(m)
DSGE.update!(dsgevar, observables =
collect(keys(m.observables)),
shocks = collect(keys(m.exogenous_shocks)),
lags = 4, λ = dsgevar_λ)
```

Or we could have also run

```
dsgevar = DSGEVAR(m, "ss3")
DSGE.update!(dsgevar, shocks = collect(keys(m.exogenous_shocks)), observables = collect(keys(m.observables)))
```
=#

# Estimate the DSGEVAR
if estimate_dsgevar
    if do_parallel
        my_procs = addprocs(n_workers)
        @everywhere using DSGE, OrderedCollections
    end

    estimate(dsgevar, data; run_csminwel = false, verbose = :none)

    if do_parallel
        rmprocs(my_procs)
    end
end

# Compute VAR coefficients and innovations variance-covariance matrices
# as well as population moments and the underlying DSGE's state space representation
if get_VAR_system
    var_system = Dict()
    var_system[:var_approx] = Dict()
    var_system[:dsgevar_λ] = Dict()
    DSGE.update!(dsgevar, modal_paras)

    # Coefficients and innovations variance-covariance matrices when approximating the DSGE
    # and when using the DSGE as a prior for the VAR
    var_system[:var_approx][:β], var_system[:var_approx][:Σ] = compute_system(dsgevar; use_intercept = true)
    var_system[:dsgevar_λ][:β], var_system[:dsgevar_λ][:Σ]   = compute_system(dsgevar, data)

    # Population moments
    var_system[:var_approx][:yyyyd], var_system[:var_approx][:xxyyd], var_system[:var_approx][:xxxxd] =
        compute_system(dsgevar; use_intercept = true, get_population_moments = true)
    var_system[:dsgevar_λ][:yyyyd], var_system[:dsgevar_λ][:xxyyd], var_system[:dsgevar_λ][:xxxxd] =
        compute_system(dsgevar, data; get_population_moments = true)

    # DSGE state space representation using the VAR's desired observables
    var_system[:var_approx][:system] = compute_system(dsgevar; use_intercept = true, get_system = true)
    var_system[:dsgevar_λ][:system]  = compute_system(dsgevar, data; get_system = true)
end

# Compute impulse responses
if do_modal_irf
    modal_irfs = Dict()

    if use_estim_output
        DSGE.update!(dsgevar, load_draws(dsgevar, :mode))
        DSGE.update!(m,       load_draws(dsgevar, :mode))
    end

    ## DSGE impulse responses
    modal_irfs[:dsge] = Dict()
    dsge_shocks = vcat(1., zeros(n_observables(m) - 1))
    system = compute_system(m)
    modal_irfs[:dsge][:cholesky] = Dict()
    modal_irfs[:dsge][:cholesky][:states], modal_irfs[:dsge][:cholesky][:obs], modal_irfs[:dsge][:cholesky][:pseudo] =
        impulse_responses(system, impulse_response_horizons(m),
                          Matrix{Float64}(I, n_observables(m), n_observables(m)),
                          dsge_shocks; restriction = :short_run, flip_shocks = false, get_shocks = false)

    # This commented code triggers a PositiveDefiniteException, likely because with the particular combination
    # of shocks and observables chosen, the long-run identification will not work.
    # modal_irfs[:dsge][:choleskyLR] = Dict()
    # modal_irfs[:dsge][:choleskyLR][:states], modal_irfs[:dsge][:choleskyLR][:obs], modal_irfs[:dsge][:choleskyLR][:pseudo] =
    #     impulse_responses(system, impulse_response_horizons(m),
    #                       Matrix{Float64}(I, n_observables(m), n_observables(m)),
    #                       dsge_shocks; restriction = :long_run, flip_shocks = false, get_shocks = false)

    modal_irfs[:dsge][:maxBC] = Dict()
    modal_irfs[:dsge][:maxBC][:states], modal_irfs[:dsge][:maxBC][:obs], modal_irfs[:dsge][:maxBC][:pseudo] =
        impulse_responses(system, impulse_response_horizons(m), (2π / 32, 2π / 6), 1;
                          flip_shocks = false, get_shocks = false)

    ## DSGEVAR impulse responses (λ = dsgevar_λ and λ = ∞)
    modal_irfs[:var_approx] = Dict() # λ = ∞
    modal_irfs[:dsgevar_λ] = Dict()  # λ = dsgevar_λ
    for method in [:cholesky, :choleskyLR, :maxBC]
        modal_irfs[:var_approx][method] =
            impulse_responses(dsgevar, method, 1; horizon = impulse_response_horizons(dsgevar),
                              flip_shocks = false, use_intercept = true)
        modal_irfs[:dsgevar_λ][method] =
            impulse_responses(dsgevar, data, method, 1; horizon = impulse_response_horizons(dsgevar),
                              flip_shocks = false)
    end
    DSGE.update!(dsgevar; λ = Inf)
    modal_irfs[:var_approx][:rotation] =
        impulse_responses(dsgevar, data; horizon = impulse_response_horizons(dsgevar),
                          deviations = true, flip_shocks = false)
    DSGE.update!(dsgevar; λ = dsgevar_λ)

    if compare_modal_irf
        DSGE.update!(dsgevar; λ = Inf)
        modal_irfs[:var_approx][:rotation] =
            impulse_responses(dsgevar, data; horizon = impulse_response_horizons(dsgevar),
                              deviations = true, flip_shocks = false, normalize_rotation = true)
        DSGE.update!(dsgevar; λ = dsgevar_λ)
        dsge_system = compute_system(m)
        dsge_states, dsge_obs, dsge_pseudo = impulse_responses(m, dsge_system)
        plots_dict = Dict{Symbol, Dict{Symbol,Plots.Plot}}()
        for (shock, j) in m.exogenous_shocks
            plots_dict[shock] = Dict{Symbol, Plots.Plot}()
            for (obs, i) in DSGE.get_observables(dsgevar)
                plots_dict[shock][obs] = plot(1:impulse_response_horizons(m), dsge_obs[i, :, j],
                                              label = "DSGE", color = :red, linestyle = :dash, linewidth = 3)
                plot!(1:impulse_response_horizons(m), modal_irfs[:var_approx][:rotation][i, :, j],
                      label = "DSGE-VAR (\\lambda = Inf)", color = :black,
                      left_margin = 20px, linewidth = 3)
            end
        end
    end
end

if do_full_band_irf
    if do_parallel
        my_procs = addprocs(n_workers)
        @everywhere using DSGE, OrderedCollections
    end

    full_band_irfs = Dict()
    modal_irfs_wrapper = Dict() # this does the same thing as the code block above, but here, we use the wrapper script to do it.

    paras = if use_estim_output
        load_draws(dsgevar, :full)
    else
        # This code repeat the modal parameters
        # and reshapes to a ndraws x nparameters Matrix.
        convert(Matrix{Float64}, repeat(map(θ -> θ.value, DSGE.get_parameters(dsgevar)), 1, 10)')
    end
    if use_estim_output
        DSGE.update!(dsgevar, load_draws(dsgevar, :mode))
    end

    ## Choose the index of the observable corresponding to the desired orthogonal shock
    n_obs_shock = 1

    ## DSGE impulse responses
    full_band_irfs[:dsge] = Dict()
    modal_irfs_wrapper[:dsge] = Dict()
    for method in [:cholesky, :maxBC]
        full_band_irfs[:dsge][method] =
            impulse_responses(m, paras, :full, method, n_obs_shock; parallel = do_parallel,
                              flip_shocks = false, create_meansbands = create_meansbands,
                              forecast_string = forecast_string)
        modal_irfs_wrapper[:dsge][method] =
            impulse_responses(m, map(θ -> θ.value, DSGE.get_parameters(dsgevar)), :mode, method, n_obs_shock; parallel = false,
                              flip_shocks = false, create_meansbands = false,
                              forecast_string = forecast_string)
    end

    ## DSGEVAR impulse responses (λ = dsgevar_λ and λ = ∞)
    full_band_irfs[:var_approx] = Dict()
    full_band_irfs[:dsgevar_λ] = Dict()
    modal_irfs_wrapper[:var_approx] = Dict()
    modal_irfs_wrapper[:dsgevar_λ] = Dict()
    for method in [:cholesky, :choleskyLR, :maxBC]
        full_band_irfs[:var_approx][method] =
            impulse_responses(dsgevar, paras, :full, method, n_obs_shock; parallel = do_parallel,
                              use_intercept = true, flip_shocks = false, create_meansbands = create_meansbands,
                              forecast_string = forecast_string)

        # The commented code below does the same as the previous lines in the case that you do
        # not want to construct a DSGEVAR object but still want to approximate a DSGE with a VAR.
        # See the docstring for more information about the differences between the two functions.
        # full_band_irfs[:var_approx][method] =
        #     impulse_responses(m, paras, :full, method, 4, [:obs_gdp, :obs_cpi, :obs_nominalrate],
        #                       collect(keys(m.exogenous_shocks)), n_obs_shock; parallel = do_parallel,
        #                       use_intercept = true, flip_shocks = false, create_meansbands = create_meansbands,
        #                       forecast_string = forecast_string)

        modal_irfs_wrapper[:var_approx][method] =
            impulse_responses(dsgevar, paras[1, :], :mode, method, n_obs_shock; parallel = do_parallel,
                              use_intercept = true, flip_shocks = false, create_meansbands = create_meansbands,
                              forecast_string = forecast_string)

        full_band_irfs[:dsgevar_λ][method] =
            impulse_responses(dsgevar, paras, data, :full, method; parallel = do_parallel,
                              flip_shocks = false, create_meansbands = create_meansbands,
                              forecast_string = forecast_string)

        modal_irfs_wrapper[:dsgevar_λ][method] =
            impulse_responses(dsgevar, paras[1, :], data, :mode, method; parallel = do_parallel,
                              flip_shocks = false, create_meansbands = create_meansbands,
                              forecast_string = forecast_string)
    end
    DSGE.update!(dsgevar; λ = Inf)
    full_band_irfs[:var_approx][:rotation] =
        impulse_responses(dsgevar, paras, data, :full, :rotation; parallel = do_parallel,
                          flip_shocks = false, create_meansbands = create_meansbands,
                          deviations = true,
                          forecast_string = forecast_string)
    DSGE.update!(dsgevar; λ = dsgevar_λ)

    DSGE.update!(dsgevar; λ = Inf)
    modal_irfs_wrapper[:var_approx][:rotation] =
        impulse_responses(dsgevar, paras[1, :], data, :mode, :rotation; parallel = do_parallel,
                          flip_shocks = false, create_meansbands = create_meansbands,
                          deviations = true,
                          forecast_string = forecast_string)
    DSGE.update!(dsgevar; λ = dsgevar_λ)

    if do_parallel
        rmprocs(my_procs)
    end
end
