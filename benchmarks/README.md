# Benchmarking Suite

This directory contains a set of performance benchmarks for the major components of the
DSGE.jl package. It is used internally to ensure that any changes made to one
component of the codebase do not detract from the performance of any other components.

Because the performance speedwise (in terms of `times`) will depend on the machine that the
benchmarks are run on, be sure to use `overwrite_ref_trial` with a fresh trial of the
function you are benchmarking before making any changes to the source code such that your
benchmarks are compared on the same hardware specifications.

Also, the assumed directory structure is that each of your benchmarked files live within
their own directories in the top-level benchmark directory, where the reference directory
also lives. This is hard coded for now, but may change in the future.
