# DSGE.jl v0.3.2 Release Notes

## New features

- Added plotting functions: `plot_parameters`, `plot_history_and_forecast`,
  `plot_forecast_comparison`, `hair_plot`, `plot_shock_decomposition`, and
  `plot_impulse_response`
- Added `nelder_mead` optimizer

## Breaking changes

- Changed input data file names: see `get_data_filename`
  + Added dataset identifier `Setting` with key `data_id`
  + Changed `cond_id` from `Setting{String}` to `Setting{Int}`
  + Moved raw input data files from `inpath(m, "data")` to `inpath(m, "raw")`
- Added `:marginal_L` (marginal likelihood) field to `Kalman` type
- Pluralized forecast output classes `:states`, `:shocks`, and `:stdshocks`


# DSGE.jl v0.3.1 Release Notes

## Bug fixes

- Added the following subspecs:
  + Model 990, subspec 3: fixes bugs 1-4 in
    [FRBNY-DSGE/DSGE-2015-Apr#1](https://github.com/FRBNY-DSGE/DSGE-2015-Apr/issues/1)
  + Model 1002, subspec 10: corrects the definition of `betabar` to use
    `m[:σ_c]` instead of `σ_ω_star`
  + Model 1010, subspec 20: similarly corrects the definition of `betabar`

## Deprecation fixes

- Implemented `transpose` for `Parameter`s so that matrix division (i.e. the
  `(\)` operator) no longer throws a warning


# DSGE.jl v0.3.0 Release Notes

## New features

- `detexify` function turns unicode characters into ASCII strings
  before writing them to CSV.

## Breaking changes

- Changed `Dict`s of indices in model object to `OrderedDict`s
- Upgrade all code for use with Julia v0.5.1 or higher


# DSGE.jl v0.2.0 Release Notes

## New features

- Added the An-Schorfheide model, a simple three-equation New Keynesian model.
- Added Model 1010, an updated version of Model 1002.
- Added three optimization methods: `:simulated_annealing`, `:LBFGS`, and
  `:combined_optimizer` (which alternates between simulated annealing and
  LBFGS).
- Added the `PseudoObservable` type and the `pseudo_measurement` function, which
  defines pseudo-observables (linear combinations of states which are not
  observed) for each model, e.g. the output gap.
- Implemented the forecast step, a suite of functions that forecast using
  estimated parameters and compute means and bands of the forecasted series. The
  top-level functions are `forecast_one` and `means_bands_all`; see the
  [forecasting](http://frbny-dsge.github.io/DSGE.jl/latest/forecast.html) and
  [means and bands](http://frbny-dsge.github.io/DSGE.jl/latest/means_bands.html)
  for more details.

## Breaking changes

- Added the `Observable` type; replaced the `data_series` and `data_transforms`
  fields in the model type definitions with
  `observable_mappings::OrderedDict{Symbol, Observable}`, which is initialized
  in `init_observable_mappings!`.
- `kalman_filter` has been broken out into
  [StateSpaceRoutines.jl](https://github.com/FRBNY-DSGE/StateSpaceRoutines.jl).
- `estimate` now saves only parameter draws, not the associated state-space
  matrices or the last filtered states for each draw.


# DSGE.jl v0.1.5 Release Notes

## New Features

- Added Model 1002, an updated version of Model 990.
- Added documentation for Model 1002 at docs/DSGE_Model_Documentation_1002.pdf.
  This pdf includes an overview of the economic theory underlying the model, a
  summary of the model's main equations, a description of
  the data used, a table of priors for the model's parameters,
  and more.

## Deprecation Fixes

- Optim.jl's `MultivariateOptimizationResults` type requires `f_increased` field
- `MersenneTwister` must be constructed with a seed


# DSGE.jl v0.1.4 Release Notes

## New Features

- Gensys no longer throws an error when system is indeterminate;
  instead, a warning is printed to the screen.

## Bug Fixes

- Fix `OptimizationTrace` constructor according to Optim v0.6. See #6.


# DSGE.jl v0.1.3 Release Notes

## New features

- Automatic dataset download and generation
- More robust and flexible treatment of dataset- and model-related dates
- Refactored computational settings
- Improved infrastructure for organizing input/output files
- Bug fix in treatment of zero lower bound in posterior computation
- Improved test coverage and documentation

## Breaking changes

- Input data matrices are CSV instead of HDF5
- Estimation output matrices are *not* flattened when saved
