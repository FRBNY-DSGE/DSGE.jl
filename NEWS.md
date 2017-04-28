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
