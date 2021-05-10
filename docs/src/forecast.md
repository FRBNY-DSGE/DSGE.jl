# [Forecasting](@id forecast-step)

This page describes forecasting without regime-switching.
Click [here](@ref regime-switch-forecast) to learn about forecasting
with regime-switching.

```@meta
CurrentModule = DSGE
```

## Procedure

In the forecast step, we compute smoothed histories, forecast, compute shock
decompositions, and compute impulse response functions (IRFs) for states,
observables, shocks, and pseudo-observables. To run a forecast on one
combination of input parameter type (e.g. modal parameters or full-distribution)
and conditional type, call `forecast_one`. The forecast output
is written to the saveroot specified by the model object and can be loaded
with the function `read_forecast_output`.

**Main Steps:**

- *Prepare forecast inputs:* Add required output types, load data, and load
  draws of parameter vectors saved from the estimation step.

- *Compute forecast outputs:* Carry out desired combination of smoothing,
  forecasting, computing shock decompositions, and computing IRFs. See
  [Forecast Outputs](@ref) for a list of possible forecast
  outputs.

- *Save forecast outputs:* Save each forecast output as an array to its own
  file, along with some metadata.

```@docs
DSGE.forecast_one
```

For example, to do an unconditional forecast of states and observables using the
modal parameters, call:

```julia
m = AnSchorfheide()
forecast_one(m, :mode, :none, [:forecaststates, forecastobs])
```

**Full-Distribution Forecasts:**

Full-distribution forecasts are computed in blocks. The size of each block
defaults to 5000 draws (before thinning by `get_setting(m, :forecast_jstep)`),
but can be set using the `:forecast_block_size` `Setting`. For each block, draws
are read in on the originator process, then computation proceeds in parallel
using `pmap`. When all draws in the block are finished, the forecast outputs are
reassembled on the originator process and appended to the HDF5 dataset in their
respective output files.

To fully take advantage of the parallelization, the user is responsible for
adding processes before calling `forecast_one`, either by calling `addprocs` or
using one of the functions defined in
[ClusterManagers.jl](https://github.com/JuliaParallel/ClusterManagers.jl). For
example, to run a full-distribution unconditional forecast using 10 processes:

```julia
my_procs = addprocs(10)
@everywhere using DSGE

m = AnSchorfheide()
forecast_one(m, :full, :none, [:forecaststates, forecastobs])

rmprocs(my_procs)
```

Notice that it is necessary to load DSGE on all processes using
`@everywhere using DSGE` before calling `forecast_one`. It is
also sometimes necessary to load OrderedCollections on all processes
using `@everywhere using OrderedCollections`.

By default, full-distribution forecasts start from the first block. However, if
you want to start the forecast from a later block, you can also do so. For
example:

``` julia
m <= Setting(:forecast_start_block, 2,
    "Block at which to resume forecasting (possibly null)")
```

## Forecast Outputs

A forecast output (i.e. an `output_var`) is a combination of what we call a
"product" and a "class". The possible classes are states (`:states`), observables
(`:obs`), pseudo-observables (`:pseudo`), and standardized (`:stdshocks`) and
unstandardized shocks (`:shocks`). The possible forecast products are:

- Smoothed histories (`:hist`): use the smoother specified by
  `forecast_smoother(m)` to get smoothed histories of each class.

- Forecasts (`:forecast`): iterate the state space forward from the last
  filtered state, either using a specified set of shock innovations or by
  drawing these from a distribution. Forecasts in which we enforce the zero
  lower bound are denoted as `:bddforecast`.

- Shock decompositions (`:shockdec`): decompose a forecast into the
  contributions from each individual shock (e.g. TFP or monetary policy) or
  from groups of shocks by exploiting the linearity of the state space system.
  Starting from an initial state of zero,
  iterate the state space forward from the first historical period up through
  the last forecast horizon. To infer the contributions of shocks,
  use the smoothed historical shocks for one shock at a time during the historical
  periods. If a modal shock decomposition is computed, then we assume no shocks
  occur during the forecast periods. If a full-distribution shock decomposition
  is computed, then we allow shocks to occur during the forecast periods.

- Deterministic trends (`:dettrend`): iterate the state space forward from the first
  historical state up through the last forecast horizon without any shocks.

- Trends (`:trend`): for each class, just the constant term in that class's
  equation, i.e. the `CCC` vector from the transition equation for states, the
  `DD` vector from the measurement equation for observables, and the `DD_pseudo`
  vector from the pseuodo-measurement equation for pseudo-observables.

- IRFs (`:irf`): see
  [Impulse response](https://en.wikipedia.org/wiki/Impulse_response). Our IRFs are
  in response to a shock of size -1 standard deviation.

An `output_var` is then just a `Symbol` with a product and class concatenated,
e.g. `:histstates` for smoothed historical states.

It is not necessary to compute all forecast outputs in one call to
`forecast_one`. Which steps are run depends on which `output_vars` are passed
in. Most users will just want to calculate `:forecastobs`, `:histobs`,
`:forecastpseudo`, and `:histpseudo`.


### [Calculating Shock Decompositions](@id calc-shock-dec)
The user should note that shock decompositions **are not** the same thing
as forecast decompositions in DSGE.jl. The former decomposes a forecast
into the contributions of individual or groups of shocks. The latter
decomposes the change between two different forecasts into "data revision",
"news" and "re-estimation" components.

If the user wants to plot shock decompositions to decompose
a forecast's deviations from trend into the contributions of individual or groups of shocks,
then the user needs to calculate the `shockdec`, `detttrend`, and `trend` output variables.
Otherwise, the plotting script `plot_shock_decompositions` will fail.

Finally, note that full-distribution shock decompositions are very memory intensive. For `Model1002`,
we typically need 30GB-40GB using the default settings. To avoid memory problems,
the user should decrease the setting `:forecast_block_size` to a smaller number than
the default 5000. A smaller block size means that fewer draws from the posterior
are loaded into memory at any one time, hence the size of the matrices
associated with shock decompositions will also be smaller.

## Preparing Forecast Inputs

**Adding Required `output_var`s:**

The user specifies what output to compute by passing a `Vector{Symbol}` of
correctly named output variables into `forecast_one`. Underneath the hood,
we will add any additional output variables that are required
to complete a forecast using `add_requisite_output_vars`:

- If `:forecast<class>` is in `output_vars`, then `:bddforecast<class>` is also
  added. Hence we always forecast both with and without enforcing the ZLB.
- If `:shockdec<class>` is in `output_vars`, then `:dettrend<class>` and
  `:trend<class>` are also added. This is because to plot shock decompositions,
  we also need the trend and the deterministic trend.

**Loading Data:**

This is done the usual way, using `load_data` with the appropriate `cond_type`.

Note that if you are running a conditional forecast, then you should check that
the appropriate observables are found in `cond_full_names(m)` and/or
`cond_semi_names(m)`. If the observables on which you are conditioning do not
match the ones found in these settings exactly, then conditional forecasting
will not generate the results you want. To change these observables to, for example,
condition on `custom_obs1` and `custom_obs2`, add the line

```
m <= Setting(:cond_full_names, [custom_obs1, custom_obs2])
```

to your script for full conditional forecasts.
The syntax is similar for semi-conditional forecasts.
If you have already created a conditional dataset successfully, then
you should be sure to re-create the data set by running

```
load_data(m; try_disk = false, cond_type = cond_type)
```

If you do not have `try_disk = false`, then the code may load
a saved data set that does not have the conditional observables you want to use.


**Loading Draws:**

By default, the draws are loaded from the file whose path is given by
`get_forecast_input_file`. However, you can directly pass a matrix
of posterior parameter draws or a vector of the modal parameters
via the keyword argument `params` for `forecast_one`,
or you can override the default input file for
a given input type by adding entries to the `Dict{Symbol, ASCIIString}` returned
from `forecast_input_file_overrides(m)`. For example:

```julia
overrides = forecast_input_file_overrides(m)
overrides[:mode] = "path/to/input/file.h5"
```

Note that `load_draws` expects an HDF5 dataset called either `params` (for
`input_type in [:mode, :mean]`) or `mhparams` (for `input_type in [:full, :subset]`)
if the sampling method is `:MH`. If the sampling method is `:SMC`, then
`load_draws` expects a jld2 file which has a variable named `cloud`.

## Computing Forecast Outputs

For each draw of parameters, the forecast calculations are run by
the lower-level function `forecast_one_draw`.

This function is also useful on its own
when the user wants to run a single-draw forecast without writing the output to data.
A common use case is running experiments with different forecast specifications and/or alternative policy rules.
Given the appropriate inputs, `forecast_one_draw` will return a `Dict` whose keys are the names of
output variables (e.g. `:forecastobs`) and values are the corresponding matrices.
This function does not perform transformations, so the units of the output are all model units,
which are typically at the quarterly frequency, such as quarterly per-capita GDP growth.


The computations that `forecast_one_draw` can run are:

**Smoothing:**

Smoothing is necessary if either:

- You explicitly want the smoothed histories, or
- You want to compute shock decompositions or deterministic trends, which use
  the smoothed historical shocks

It is not necessary to keep track of these cases, however - `forecast_one` will
deduce from the specified `output_vars` whether or not it is necessary to filter
and smooth in order to produce your `output_vars`.

**Forecasting:**

Forecasting begins from the last filtered historical state, which is obtained
from the Kalman filter. `forecast` accepts a keyword argument `enforce_zlb`,
which indicates whether to enforce the zero lower bound. If `enforce_zlb =
true`, then if in a given period, the forecasted interest rate goes below
`forecast_zlb_value(m)`, we solve for the interest rate shock necessary to push
it up to the ZLB. A forecast in which the ZLB is enforced corresponds to the
product `:bddforecast`.

**Shock Decompositions, Deterministic Trends, and Trends:**

Since shock decompositions have an additional dimension (e.g. `nstates` x
`nperiods` x `nshocks` for a single draw of state shock decompositions, compared
to `nstates` x `nperiods` for a single draw of forecasted states), we usually
wish to truncate some periods before returning. This behavior is governed by the
`Settings` `:shockdec_starttdate` and `:shockdec_enddate`, which are of type
`Nullable{Date}`.

Deterministic trends are also saved only for `date_shockdec_start(m)` and
`date_shockdec_end(m)`. Trends are not time-dependent.

Note that shock decompositions are memory intensive. For a typical forecast
of `Model1002` without shock decompositions, only 1-2GB of memory are needed.
For a forecast with shock decompositions, roughtly 14-16GB will be needed.

**Impulse Response Functions:**

Like shock decompositions, IRFs have three dimensions (e.g. `nstates` x
`nperiods` x `nshocks`) for each draw.


## Saving Forecast Outputs

Forecast outputs are saved in the location specified by
`get_forecast_output_files(m)`, which is typically a subdirectory of
`saveroot(m)`. Each `output_var` is saved in its own JLD file, which contains
the following datasets:

- `arr::Array`: actual array of forecast outputs. For trends, this array is of
  size `ndraws` x `nvars`. For histories, forecasts, and deterministic trends,
  it is `ndraws` x `nvars` x `nperiods`. For shock decompositions and IRFs, it
  is `ndraws` x `nvars` x `nperiods` x `nshocks`. (In all of these, `nvars`
  refers to the number of variables of the output class.)

- `date_indices::Dict{Date, Int}`: maps `Date`s to their indices along the
  `nperiods` dimension of `arr`. Not saved for IRFs.

- `<class>_names::Dict{Symbol, Int}`: maps names of variables of the output
  class (e.g. `:OutputGap`) into their indices along the `nvars` dimension of
  `arr`.

- `<class>_revtransforms::Dict{Symbol, Symbol}`: maps names of variables to the
  names of the reverse transforms (from model units into plotting units)
  associated with those variables. For example,
  `pseudoobservable_revtransforms[:π_t] = :quartertoannual`.

- `shock_names::Dict{Symbol, Int}`: for shock decompositions and IRFs only, maps
  names of shocks into their indices along the `nshocks` dimension of `arr`.

Some helpful functions for getting file names, as well as reading and writing
forecast outputs, include:

- `get_forecast_input_file`
- `get_forecast_filename`
- `get_forecast_output_files`
- `write_forecast_outputs`
- `write_forecast_block`
- `write_forecast_metadata`
- `read_forecast_metadata`
- `read_forecast_output`

## Forecasting Functions
```@autodocs
Modules = [DSGE]
Pages   = ["drivers.jl", "forecast.jl", "smooth.jl", "shock_decompositions.jl", "io.jl"]
Order   = [:function]
```
