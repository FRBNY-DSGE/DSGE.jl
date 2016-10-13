# Running the Code

## Running with Default Settings

To run the estimation step in Julia, simply create an instance of the model object and pass
it to the `estimate` function -- see an [example](doc/examples/run_default.jl).

```julia
# construct a model object
m = Model990()

# estimate as of 2015-Q3 using the default data vintage from 2015 Nov 27
m <= Setting(:data_vintage, "151127")
m <= Setting(:date_mainsample_end, quartertodate("2015-Q3"))

# reoptimize parameter vector, compute Hessian at mode, and full posterior
# parameter sampling
estimate(m)

# produce LaTeX tables of parameter moments
compute_moments(m)
```

By default, the `estimate` routine loads the dataset, reoptimizes the initial parameter
vector, computes the Hessian at the mode, and conducts full posterior parameter sampling.
(The initial parameter vector used is specified in the model's constructor.)

To use updated data or alternative user-specified datasets, see [Input Data](#input-data).

The user may want to avoid reoptimizing the parameter vector and calculating the
Hessian matrix at this new vector. Please see [Reoptimizing](#reoptimizing)
below.

For more details on changing the model's default settings, parameters, equilibrium
conditions, etc., see [Advanced usage](#advanced-usage).

## Input/Output Directory Structure

The *DSGE.jl* estimation uses data files as input and produces large data files
as outputs. One estimation saves several GB of parameter draws and
related outputs. It is useful to understand how these files are loaded/saved
and how to control this behavior.

### Directory Tree
The following subdirectory tree indicates the default locations of
these input and outputs. Square brackets indicate directories in the tree that
will become relevant as future features are implemented.
- `<dataroot>/`: Root data directory.
  - `data/`:  Macroeconomic input data series.
  - `cond/`: Conditional data, i.e.
    ["nowcast"](https://en.wikipedia.org/wiki/Nowcasting_%28economics%29).
  - `user/`: User-created or sample model input files.

- `<saveroot>/`: Root save directory.
  - `output_data/`
    - `m990/`: Input/output files for the `Model990` type. A model of type
      `SPEC` will create its own save directory `SPEC/` at this  level in the
      directory tree.
      - `ss0/`: Subdirectory for subspec 0. A model of a different subspec will have similar
          directories at this level of the tree.
        - `estimate/`
          - `figures/`: Plots and other figures
          - `tables/`: LaTeX tables
          - `raw/`: Raw output data from estimation step
          - `work/`: Derived data files created using `raw/` files as input
        - [`xxx/`]: Other model outputs, such as forecasts, impulse response
          functions, and shock decompositions; subdirectory structure mirrors that of
          `estimate`.

### Directory Paths

By default, input/output directories are located in the *DSGE.jl* package, along
with the source code. Default values of the input/output directory roots:
- `saveroot(m)`: `"$(Pkg.dir())/DSGE/save"`
- `dataroot(m)`: `"$(Pkg.dir())/DSGE/save/input_data"`

Note these locations can be overridden as desired. See [Settings](#model-settings) below for more
details.
```julia
m <= Setting(:saveroot, "path/to/my/save/root")
m <= Setting(:dataroot, "path/to/my/data/root")
```

Utility functions are provided to create paths to input/output files. These should be used
for best results.
- `inpath`: Return path to directory/file for *input data*, input conditional data, and
    user-provided sample files, respectively. See `?inpath` for more details.
    - Try `inpath(m, "data")`, a helpful call that displays the containing directory of
      input data files.
- `rawpath`: Return path to directory/file for a given output type for *raw model output*.
    See `?rawpath` for more details.
- `workpath`: Return path to directory/file for a given output type for *transformed* or
    *intermediate model output*. See `?workpath` for more details.
- `tablespath`: Return path to directory/file for a given output type for *results tables*
    or other *textual results*.
- `figurespath`: Return path to directory/file for a given output type for *results
    figures*.

