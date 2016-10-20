# Forecasting

```@meta
CurrentModule = DSGE
```

## Procedure

In the forecast step, we compute smoothed histories, forecast, and compute shock
decompositions for states, observables, shocks, and so-called
pseudo-observables. (Like observables, pseudo-observables are linear
combinations of the states, but we may not necessarily have historical data for
them. The canonical example is the output gap, the difference between GDP and
potential output.) The main function used in the forecast step is `forecast_one`.

**Main Steps:**

- *Prepare forecast inputs:* Load draws of parameter vectors, transition
  matrices, and states from Metropolis-Hastings and reshape them for use in the
  next steps.

- *Filter and smooth histories:* Filter and smooth to get smoothed historical
  states and shocks.

- *Forecast:* Iterate the state space forward from the last historical state.

- *Compute shock decompositions:* Iterate the state space forward from the first
  historical state, shutting down all but one shock at a time.

It is not necessary to do all of filtering and smoothing, forecasting, and
computing shock decompositions in one call to `forecast_one`. Which steps are
run depends on the keyword argument `output_vars` passed in. To forecast using
one combination of input type, conditional data type, and output variables,
specify the keyword arguments as below:

```@docs
DSGE.forecast_one
```

For example, to do an unconditional forecast of states and observables using the
modal parameters, call:

```julia
m = Model990()
df = load_data(m)
forecast_outputs = forecast_one(m, df; input_type = :mode, cond_type = :none,
                       output_vars = [:forecaststates, forecastobs])
```

**Full-Distribution Forecasts**

Since Julia is just-in-time compiled, a function `f` is compiled only on the
first call to `f` in a given Julia session. On this first call, `f` both takes
longer to run and allocates more memory than in future calls. (See
[Performance Tips](http://docs.julialang.org/en/latest/manual/performance-tips/)
in the Julia manual for more details.) For this reason, it is important to
compile `forecast_one` by running it once with `input_type = :full` on a small
set of draws before running it again on your full set of draws. The function
[`compile_forecast_one`](@ref) is provided to simplify this step; call it with
the same keyword arguments (except omitting `input_type`, because
`input_type = :full` is assumed) that you will call `forecast_one`. For example:

```julia
m = Model990()
df = load_data(m)
compile_forecast_one(m, df; cond_type = :none, output_vars = [:forecaststates, :forecastobs])
forecast_outputs = forecast_one(m, df; input_type = :full, cond_type = :none,
    output_vars = [:forecaststates, :forecastobs])
```


## Preparing Forecast Inputs

The function `prepare_forecast_inputs` reads in the draw or draws of parameters,
transition matrices, and last-historical-period filtered states (all saved by
Metropolis-Hastings during the estimation step). The saved transition matrix
draws are reshaped from the format in which they were saved and assembled into
`System` objects.

```@docs
DSGE.prepare_forecast_inputs
```

By default, the draws are loaded from the file whose path is given by
`get_input_file(m, input_type)`. However, you can override the default input
file for a given input type by adding entries to the `Dict{Symbol, ASCIIString}`
returned from `forecast_input_file_overrides(m)`. For example:

```julia
overrides = forecast_input_file_overrides(m)
overrides[:mode] = "path/to/input/file.h5"
```

`prepare_forecast_inputs` calls `load_draws`, `prepare_systems`, and
`prepare_states`:

```@docs
DSGE.load_draws
DSGE.prepare_systems
DSGE.prepare_states
```


## Filtering and Smoothing

In `filterandsmooth_all`, we use the Kalman filter on all draws to get filtered
historical states and shocks, and then use the smoother specified by
`forecast_smoother(m)` to get smoothed states and shocks. We then apply the
measurement and pseudo-measurement equations to get smoothed observables and
pseudo-observables, respectively.

Filtering and smoothing is necessary in certain situations:

- If you explicitly want the smoothed histories
- If you want to compute shock decompositions, which use the smoothed historical
  shocks
- If you want to forecast using conditional data. This is because the
  "historical data" loaded by `load_data` when
  `cond_type in [:semi, :full]` actually includes partial data for
  `date_forecast_start(m)` and later. The smoothed histories in these so-called
  conditional periods are therefore actually forecasts. Since we start our
  forecasts by iterating forward from the last "historical" state, we must
  transplant the conditional periods from the smoothed histories to the
  forecasts, so that in the end the smoothed histories end in
  `date_mainsample_end(m)` and the forecasts begin at
  `date_forecast_start(m)`.

It is not necessary to keep track of these cases, however - `forecast_one` will
deduce from the specified `output_vars` whether or not it is necessary to filter
and smooth in order to produce your `output_vars`.

To filter and smooth over many draws, call `filterandsmooth_all`, which in turn
calls `filterandsmooth` on each draw. You can also just call `filterandsmooth`
on one draw.

```@docs
DSGE.filterandsmooth_all
DSGE.filterandsmooth
```


## Forecasting

In `forecast`, we iterate the state space forward from the last
filtered state, either using a specified set of shock innovations or by
drawing these from a distribution. The end filtered states for each draw were
obtained in [`prepare_forecast_inputs`](@ref) in one of two ways:

- If `cond_type = :none`, [`load_draws`](@ref) reads in the end filtered states
  saved in Metropolis-Hastings.
- If `cond_type in [:semi, :full]`, the saved end filtered states correspond to
  `date_mainsample_end(m)` (the last main-sample period), but we wish to start
  forecasting from `date_conditional_end(m)` (the last conditional data
  period). Hence [`prepare_states`](@ref) re-runs the Kalman filter in order to
  get the filtered states in `date_conditional_end(m)`.

After forecasting states by iterating forward from the last filtered state, we
again apply the measurement and pseudo-measurement equations to get the
forecasted observables and pseudo-observables.

As in the previous step, you have the option of forecasting over many draws with
`forecast` or using one draw with `compute_forecast`:

```@docs
DSGE.forecast
DSGE.compute_forecast
```


## Shock Decompositions

To get shock decompositions, we iterate the state space forward from the first
historical state up through the last forecast horizon, using the smoothed
histories for one shock at a time. As before, we decompose the observables and
pseudo-observables by apply the respective measurement equations.

Since shock decompositions have an additional dimension (e.g. `nstates` x
`nperiods` x `nshocks` for a single draw of state shock decompositions, compared
to `nstates` x `nperiods` for a single draw of forecasted states), we usually
wish to truncate some periods before returning. This behavior is governed by the
`Settings` `:shockdec_starttdate` and `:shockdec_enddate`, which are of type
`Nullable{Date}`.

Shock decompositions are computed over many draws using `shock_decompositions`,
and on one draw using `compute_shock_decompositions`:

```@docs
DSGE.shock_decompositions
DSGE.compute_shock_decompositions
```


## Distributed Arrays

Distributed arrays (`DArray`s), provided by
[DistributedArrays.jl](https://github.com/JuliaParallel/DistributedArrays.jl),
are a datastructure in which large arrays are distributed over several
processes, reducing the memory load on any single process. This is important for
full-distribution forecasts because we perform filtering and smoothing,
forecasting, and computing shock decompositions over a large number of draws,
which is very memory-intensive. Draws are distributed in
`prepare_forecast_inputs`, and remain distributed throughout the forecast step,
including while saving forecast outputs to file.

The keyword argument `procs` to `forecast_one` and other forecast-step functions
specifies the process IDs over which to distribute the various `DArrays`. The user
is responsible for adding processes before calling `forecast_one`, either by
calling `addprocs` or using one of the functions defined in
[ClusterManagers.jl](https://github.com/JuliaParallel/ClusterManagers.jl). For
example, to run a full-distribution unconditional forecast using 10 processes:

```julia
my_procs = addprocs(10)
@everywhere using DSGE

m = Model990()
df = load_data(m)
m <= Setting(:use_parallel_workers, true)
forecast_one(m, df; input_type = :full, cond_type = :none,
    output_vars = [:forecaststates, forecastobs], procs = my_procs)

rmprocs(my_procs)
```

Notice that it is necessary to load DSGE on all processes using
`@everywhere using DSGE` before calling `forecast_one`.

If `procs` is not provided, the default value is `[myid()]`, the current process
ID. Moreover, before distributing draws, `forecast_one` (or any other function
taking the keyword argument `procs`) calls [`reset_procs`](@ref), which will reset
`procs` to `[myid()]` if one of:


- `input_type` is one of `[:init, :mode, :mean]`. That is, we have only one draw
  and hence don't want to distribute it over many processes.
- `use_parallel_workers(m) = false`, i.e. if the user has not specifically
  updated the `Setting` `:use_parallel_workers` to `true`.

Furthermore, DistributedArrays requires that the dimensions of a given `DArray`
distribute evenly over the given `procs`. Since we only distribute `DArray`s
over the draw dimension, this means that the number of draws (after thinning in
the forecast step) must be divisible by the number of `procs`. If this is not
the case, a warning is thrown and the draws are truncated to ensure
divisibility. (See [`load_draws`](@ref) for further clarification.)

Each forecast output (for example, forecasts of observables) is returned from
its respective function as a `DArray`, and each of these is written to file
without ever converting back to an `Array` (that is, without copying the local
parts of the `DArray` back to the originator process). I/O with `DArrays` in
DSGE is handled by `write_darray` and `read_darray`:

```@docs
DSGE.write_darray
DSGE.read_darray
```


## Forecast Utilities

Other functions used during the forecast step include:

```@docs
DSGE.compute_system
DSGE.reset_procs
DSGE.get_jstep
DSGE.get_input_file
DSGE.get_output_vars
DSGE.get_all_output_vars
DSGE.get_output_files
DSGE.write_forecast_metadata
DSGE.compile_forecast_one
```
