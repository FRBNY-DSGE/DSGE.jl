# Forecasting

```@meta
CurrentModule = DSGE
```

## Procedure

In the forecast step, we compute smoothed histories, forecast, compute shock
decompositions, and compute impulse response functions (IRFs) for states,
observables, shocks, and pseudo-observables. To run a forecast on one
combination of input parameter type (e.g. modal parameters or full-distribution)
and conditional type, call `forecast_one`. To run several combinations, use
`forecast_all`, which iterates through all combinations, calling `forecast_one`
on each.

**Main Steps:**

- *Prepare forecast inputs:* Load draws of parameter vectors saved from the
  estimation step and solve the system for each draw. If necessary for the
  desired forecast output, load data and run the Kalman filter on each draw.

- *Smooth histories:* Get smoothed historical states and shocks.

- *Forecast:* Iterate the state space forward from the last historical state.

- *Compute shock decompositions:* Iterate the state space forward from the first
  historical state, shutting down all but one shock at a time.

- *Compute IRFs:* For each shock, calculate the response to a shock of size -1
  standard deviation.

It is not necessary to do all of filtering and smoothing, forecasting, computing
shock decompositions, and computing IRFs in one call to `forecast_one`. Which
steps are run depends on which `output_vars` are passed in.

```@docs
DSGE.forecast_all
DSGE.forecast_one
```

For example, to do an unconditional forecast of states and observables using the
modal parameters, call:

```julia
m = Model990()
forecast_outputs = forecast_one(m, :mode, :none, [:forecaststates, forecastobs])
```

To forecast states and observables at the mode for all conditional types, call:


```julia
m = Model990()
forecast_all(m, [:mode], [:none, :semi, :full], [:forecaststates, :forecastobs])
```

**Full-Distribution Forecasts**

Since Julia is just-in-time compiled, a function `f` is compiled only on the
first call to `f` in a given Julia session. On this first call, `f` both takes
longer to run and allocates more memory than in future calls. (See
[Performance Tips](http://docs.julialang.org/en/latest/manual/performance-tips/)
in the Julia manual for more details.) For this reason, it is important to
compile `forecast_one` by running it once with `input_type = :full` or
`input_type = :subset` on a small set of draws before running it again on your
full set of draws. The function [`compile_forecast_one`](@ref) is provided to
simplify this step. Call it using the same `output_vars` that you plan to call
`forecast_one` with. For example:

```julia
m = Model990()
compile_forecast_one(m, [:forecaststates, :forecastobs])
forecast_outputs = forecast_one(m, :full, :none, [:forecaststates, :forecastobs])
```


## Preparing Forecast Inputs

Before running the main functions called by `forecast_one`, we must:

1. Call `load_data` with the appropriate `cond_type`
2. Load the parameter draw or draws using `load_draws`
3. Compute state-space systems for each draw using `prepare_systems`
4. Run the Kalman filter on each draw using `filter_all` to get
   last-historical-period filtered states from which to begin the forecast

Not all of these steps are necessary for all `output_vars`; which steps get run
depends on the `output_vars` provided to `forecast_one`.

The function `prepare_forecast_inputs!` ensures that the provided keyword
arguments to `forecast_one` are well-formed, and loads the ones that are not
provided.

Why aren't all forecast inputs required arguments to `forecast_one`? For
example, if you want to run both conditional and unconditional full-distribution
forecasts, you will need to load the data and run the Kalman filter separately
for each conditional case, since conditional data has an additional
period. However, you don't need to load the draws twice, since they are not
affected by the `cond_type`. Hence you can load the draws once and pass in
`systems` as a keyword argument, but have `prepare_forecast_inputs!` call
`load_draws` and `filter_all` for you inside each `forecast_one` call:

```julia
m = Model990()
params = load_draws(m, :full)
systems = prepare_systems(m, :full, params)

for cond_type in [:none, :full]
    forecast_one(m, :full, cond_type, [:forecastobs]; systems = systems)
end
```

```@docs
DSGE.prepare_forecast_inputs!
```

By default, the draws are loaded from the file whose path is given by
[`get_input_file`](@ref). However, you can override the default input
file for a given input type by adding entries to the `Dict{Symbol, ASCIIString}`
returned from `forecast_input_file_overrides(m)`. For example:

```julia
overrides = forecast_input_file_overrides(m)
overrides[:mode] = "path/to/input/file.h5"
```

`prepare_forecast_inputs!` calls [`load_data`](@ref), `load_draws`, `prepare_systems`, and
`filter_all`:

```@docs
DSGE.load_draws
DSGE.prepare_systems
DSGE.filter_all
```


## Smoothing

In `smooth_all`, we use the smoother specified by `forecast_smoother(m)` to get
smoothed historical states and shocks. We then apply the measurement and
pseudo-measurement equations to get smoothed observables and pseudo-observables,
respectively.

Smoothing is necessary if:

- You explicitly want the smoothed histories
- You want to compute shock decompositions, which use the smoothed historical
  shocks

It is not necessary to keep track of these cases, however - `forecast_one` will
deduce from the specified `output_vars` whether or not it is necessary to filter
and smooth in order to produce your `output_vars`.

To smooth over many draws, call `smooth_all`, which in turn calls `smooth` on
each draw. You can also just call `smooth` on one draw.

```@docs
DSGE.smooth_all
DSGE.smooth
```


## Forecasting

In `forecast`, we iterate the state space forward from the last filtered state,
either using a specified set of shock innovations or by drawing these from a
distribution. For each draw, the end filtered states were obtained by running
the Kalman filter in the preparing inputs stage. After forecasting states by
iterating forward from the last filtered state, we again apply the measurement
and pseudo-measurement equations to get the forecasted observables and
pseudo-observables.

We can choose whether or not to enforce the zero lower bound by setting the
keyword argument `enforce_zlb` to `forecast` appropriately. If `enforce_zlb =
true`, then if in a given period, the forecasted interest rate goes below
`forecast_zlb_value(m)`, we solve for the interest rate shock necessary to push
it up to the ZLB. At the `forecast_one` level, we specify enforcing the ZLB by
using the `output_vars` `:forecaststatesbdd`, `:forecastobsbdd`,
`:forecastpseudobdd`, and `forecastshocksbdd`.

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


## Impulse Response Functions

See [Impulse response](https://en.wikipedia.org/wiki/Impulse_response). Our IRFs
are in response to a shock of size -1 standard deviation. Like shock
decompositions, IRFs have three dimensions (e.g. `nstates` x `nperiods` x
`nshocks`) for each draw.

As before, we can compute IRFs for all draws using `impulse_responses`, and for
one draw using `compute_impulse_response`:

``` @docs
DSGE.impulse_responses
DSGE.compute_impulse_response
```


## Distributed Arrays

Distributed arrays (`DArray`s), provided by
[DistributedArrays.jl](https://github.com/JuliaParallel/DistributedArrays.jl),
are a datastructure in which large arrays are distributed over several
processes, reducing the memory load on any single process. This is important for
full-distribution forecasts because we perform filtering and smoothing,
forecasting, computing shock decompositions, and computing IRFs over a large
number of draws, which is very memory-intensive. Draws are distributed in
`prepare_systems`, and remain distributed throughout the forecast step,
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
m <= Setting(:use_parallel_workers, true)
compile_forecast_one(m, [:forecaststates, :forecastobs]; procs = my_procs)
forecast_one(m, :full, :none, [:forecaststates, forecastobs]; procs = my_procs)

rmprocs(my_procs)
```

Notice that it is necessary to load DSGE on all processes using
`@everywhere using DSGE` before calling `forecast_one`. The same `procs` should
be used in the call to `compile_forecast_one`.

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
parts of the `DArray` back to the originator process). I/O with `DArray`s in
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
DSGE.get_forecast_output_files
DSGE.write_forecast_outputs
DSGE.write_forecast_metadata
DSGE.read_forecast_metadata
DSGE.compile_forecast_one
```
