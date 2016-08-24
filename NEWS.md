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
