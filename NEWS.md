# DSGE.jl 1.2.2
## Bug fixes
+ Fix print statements when using `forecast_one` and `verbose = :high`

# DSGE.jl 1.2.1
## Miscellaneous
+ Fix docs
+ Change default branch from master to main

# DSGE.jl 1.2.0
## Version Restrictions
+ Minimum CSV version is v0.8
+ Minimum ModelConstructors version is v0.2.4

## New features and enhancements
+ Functions for calculating k-period ahead expectations and sums for measurement equation
+ Time-varying information sets for regime-switching
+ Time-varying imperfect awareness/credibility
+ Ability to use temporary policies during history (previously only available as alternative policy during forecast)
+ Refactor alternative policy code to use regime-switching
+ Regime-switching allowed during forecast for shock decompositions, deterministic trends, and trends
+ Regime-switching estimation of DSGE models using MH and SMC

## Bug fixes and cleanup
+ Fix error in calculating the diagonal of the Hessian with parallel workers
+ Update tests to handle seed changes in Julia 1.5
+ Convert missing values to NaNs when calculating means and bands
+ Fix accidental assumption that all fixed parameters occur at the end of the parameter vector in Metropolis-Hastings
+ Fix incorrect parameter blocking in Metropolis-Hastings
+ Refactor regime-switching code to make it easier to use and maintain
+ Update syntax for HDF5 deprecations

# DSGE.jl 1.1.6
## Version Restrictions
+ Raise Julia compatibility to all v1.x.
+ Increase version restrictions for most packages to most recent ones (e.g. Plots)

## New features and enhancements
+ Implement (exogenously) regime-switching DSGE models (but estimations are not tested)
+ Implement (exogenously) regime-switching forecasts during the history and forecast horizon
+ New subspecs of `Model1002` to model the impact of the COVID-19 pandemic
+ Interface for the automated addition of anticipated shocks to DSGE models for shocks other than monetary policy shocks.
+ Interface for calculating weighted averages of different full-distribution forecasts
+ Add Nominal GDP targeting and average inflation targeting as alternative policies
+ Automatic enforcement of ZLB as a temporary alternative policy

## Bug fixes and cleanup
+ Filepathing for Windows OS should work properly now (although we still do not test on Windows OS)

# DSGE.jl 1.1.5
## New features and enhancements
+ Implement DSGE-VECMs and provide functionality to compute impulse responses.
+ Extend DSGE-VAR rotation impulse responses to compute in deviations from baseline.
+ Example for running impulse responses of a a DSGE-VECM.
+ Extend example for DSGE-VARs to allow plotting of modal impulse responses.

## Bug fixes and cleanup
+ Fix summary statistics when loading data
+ Fix bugs in the documentation.

# DSGE.jl 1.1.3
## New features and enhancements
+ Implement DSGE-VARs and provide functionality to estimate them and compute impulse responses.
+ Example for running a forecast decomposition.
+ Example for running a DSGE-VAR.

## Bug fixes and cleanup
+ Allow user to avoid running csminwel when using SMC and calling estimate.

# DSGE.jl 1.1.2
## New features and enhancements
+ Support for bridging from another estimation (in SMC.jl)--this update ensures compatibility
+ Automatically run csminwel after an SMC estimation to ensure you've found the mode with kwarg run_csminwel in smc2

# DSGE.jl 1.1.1
## New features and enhancements
+ Bug fixes and improvements
+ Add models m904 (SWFF) and m805 (SWpi)

# DSGE.jl 1.1.0
## New features and enhancements
+ Functionality to specify heterogeneous agent models, in addition to representative agent models.
+ Methods for the solution of heterogeneous agent continuous time models
+ Example model implementations

## Patches
+ Sorted out compatibility issues with `StateSpaceRoutines.jl`. Dropped compatibility with `v0.7`.

# DSGE.jl 1.0.0
## New features and enhancements
+ Package officially compatible with `v0.7`, `v1.0`, `v1.1`
+ Enable parameter blocking in Metropolis-Hastings algorithm
+ Functionality to simulate data from a model
+ Expanded test suite and example files

# DSGE.jl 0.8.1
Bug fixes and cleanup

# DSGE.jl 0.8.0

## New features and enhancements
+ Breaks out SMC and Model Constructor objects into SMC.jl and ModelConstructors.jl
+ MH can just take vector of Parameters, likelihood function, and data (like SMC but still exists only in DSGE.jl package)
+ DSGE.jl depends on SMC.jl and ModelConstructors.jl but for users who only need SMC or Model Constructor utilities, don't need to load DSGE.jl package anymore

# DSGE.jl 0.7.2

## New features and enhancements
+ Adds the ability to easily create packets of results from estimating and forecasting a DSGE model. An example script is provided in docs/examples/make_packet.jl.

# DSGE.jl 0.7.1
Bug fixes and cleanup

# DSGE.jl 0.7.0 Release Notes

## New features and enhancements
+ Adds Sequential Monte Carlo (SMC) as an alternative to Metropolis Hastings for estimating models. Latest release with all bug fixes and speed improvements

# DSGE.jl 0.6.0 Release Notes

## New features and enhancements
   + Forecast decompositions user to compare two forecasts and break down why forecasts have changed (by shock, differences in parameters, differences in data, etc.)
   + Specify size of desired impulse responses (on impact) and flip shocks
   + `filter_shocks` method allows user to obtain only filtered-only shocks (as opposed to filtered *and smoothed* shocks)
   + Compatibility with SMC (full version to be released soon)
   + New functions to deal with Date objects
   + Compatibility with likelihood-only Kalman Filter and Chandrsekhar recursions (see StateSpaceRoutines.jl) which offer large speedups

## Bug Fixes and cleanup
   + Fixes to scenarios code, inluding rectifying date labels

# DSGE.jl v0.5.1 Release Notes
+ Patch release to fix failing test

# DSGE.jl v0.5.0 Release Notes
## Breaking changes
+ Upgraded all code for use with Julia v0.7.0.
+ Changed file-saving dependency from JLD to JLD2.
+ Updated data loading and other machinery to rely on the Missing type for missing data as opposed to NaN.

# DSGE.jl v0.4.2 Release Notes
## New features and enhancements
  + Implement benchmarking suite to benchmark code performance.
  + Make DSGE compatible with improvements to the Kalman filter.

## Bug fixes and cleanup
  + Clean up unit tests.
  + Clean up travis file.
  + Fix dependency issues with NLOpt and Optim.

# DSGE.jl v0.4.1 Release Notes
## Bug fixes and cleanup
- Addressed deprecations and warnings for:
  + Non-vectorized functions (e.g. `log`, `!`)
  + `vcat` on DataFrames with different column names
  + Optimization-related functions for new releases of Optim.jl
  + Array declaration

- Tidied:
  + Plotting code
  + Meansbands computation
  + Various models
  + Scenarios code
  + Tests

# DSGE.jl v0.4.0 Release Notes

## New features

- Added `nelder_mead` optimizer
- Added forecasting under alternative policies (`AltPolicy`) and alternative
  scenarios (`AbstractScenario`)
- Added plotting functions: `plot_parameters`, `plot_history_and_forecast`,
  `plot_forecast_comparison`, `hair_plot`, `plot_shock_decomposition`,
  `plot_impulse_response`, `plot_altpolicies`, and `plot_scenario`

## Breaking changes

- Upgraded all code for use with Julia v0.6.0 or higher
- Changed input data file names: see `get_data_filename`
  + Added dataset identifier `Setting` with key `data_id`
  + Changed `cond_id` from `Setting{String}` to `Setting{Int}`
  + Moved raw input data files from `inpath(m, "data")` to `inpath(m, "raw")`
- Added `:marginal_L` (marginal likelihood) field to `Kalman` type
- Removed `MM` and `VVall` fields from `Measurement` type
- Pluralized forecast output classes `:states`, `:shocks`, and `:stdshocks`
- Stopped adding back population growth when reverse transforming shock
  decompositions and deterministic trends
- Stopped adding trends to and detrending shock decompositions and deterministic
  trends
- Changed pseudo-observable implementation to correspond one-to-one with
  observables
  + Changed `PseudoObservableMapping` type (and field in `System` type) to
    `PseudoMeasurement`
  + Added `m.pseudo_observables` and `m.pseudo_observable_mappings` fields to
    `AbstractModel` subtypes
  + Pseudo-observable-related things are no longer `Nullable`. Instead, if no
    pseudo-measurement equation is implemented, the fields in the model object
    are empty dictionaries
- Refactored means and bands computation
  + Renamed `means_bands_all` to `compute_meansbands`
  + Renamed `meansbands_matrix_all` to `meansbands_to_matrix`


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
