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
