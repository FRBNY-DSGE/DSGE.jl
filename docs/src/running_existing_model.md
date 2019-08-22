# Running an Existing Model

The DSGE.jl package provides several example models:

- A simple three-equation DSGE model from [An and Schorfheide (2006)](https://sites.sas.upenn.edu/schorf/files/er-final.pdf)
- The well-known [Smets and Wouters (2007)](https://www.aeaweb.org/articles?id=10.1257/aer.97.3.586) model
- The New York Fed DSGE model (version 990.2), which was introduced in [this blog post](http://libertystreeteconomics.newyorkfed.org/2015/05/the-frbny-dsge-model-forecast-april-2015.html)
- The New York Fed DSGE model (version 1002.9), which is documented [here](https://github.com/FRBNY-DSGE/DSGE.jl/blob/master/docs/DSGE_Model_Documentation_1002.pdf)
- The New York Fed DSGE model (version 1010.18)

You can run these models using the description provided here. If you
were to implement another model using DSGE.jl, these procedures can also be used to
estimate those models.

## Running with Default Settings

To estimate and forecast in Julia, simply create an instance of the model object
and call `estimate` and `forecast_all`. A minimal
[example](https://github.com/FRBNY-DSGE/DSGE.jl/blob/master/docs/examples/run_default.jl)
is reproduced below:

```julia
# estimate as of 2015-Q3 using the default data vintage from 2015 Nov 27
custom_settings = Dict{Symbol, Setting}(
    :data_vintage        => Setting(:data_vintage, "151127"),
    :date_forecast_start => Setting(:date_forecast_start, quartertodate("2015-Q4")))

# construct a model object
m = Model990(custom_settings = custom_settings)

# reoptimize parameter vector, compute Hessian at mode, and full posterior
# parameter sampling
estimate(m)

# produce LaTeX tables of parameter moments
compute_moments(m)

# forecast and compute means and bands using 10 processes
my_procs = addprocs(10)
@everywhere using DSGE

forecast_one(m, :full, :none, [:forecaststates, :forecastobs])
compute_meansbands(m, :full, :none, [:forecaststates, :forecastobs])
rmprocs(my_procs)
```

For more details on changing the model's default settings, parameters, equilibrium
conditions, etc., see [Advanced Usage](@ref).

By default, the `estimate` routine loads the dataset, reoptimizes the initial parameter
vector, computes the Hessian at the mode, and conducts full posterior parameter sampling
using Metropolis-Hastings.
(The initial parameter vector used is specified in the model's constructor.)
Further options for estimation are described in [Estimation](@ref estimation-step):

- To use updated data or alternative user-specified datasets, see [Input Data](@ref input-data-step).
- The user may want to avoid reoptimizing the parameter vector and calculating
  the Hessian matrix at this new vector. Please see [Reoptimizing](@ref
  estimation-reoptimizing).

For more information on the many types of forecasts that can be run on an
existing or user-defined model, see [Forecasting](@ref forecast-step).



## Input/Output Directory Structure

The *DSGE.jl* estimation uses data files as input and produces large data files
as outputs. One estimation saves several GB of parameter draws and
related outputs. It is useful to understand how these files are loaded/saved
and how to control this behavior.

### Directory Tree
The following subdirectory tree indicates the default locations of
these input and outputs. Square brackets indicate directories in the tree that
will become relevant as future features are implemented.

*Note that this directory tree is not linked, although it appears to be.*

```@contents
Pages = ["io_dirtree.md"]
Depth = 5
```

### Directory Paths

By default, input/output directories are located in the *DSGE.jl* package, along
with the source code. Default values of the input/output directory roots:

- `saveroot(m)`: `"$(Pkg.dir())/DSGE/save"`
- `dataroot(m)`: `"$(Pkg.dir())/DSGE/save/input_data"`

Note these locations can be overridden as desired. See [Model Settings](@ref) for more
details.

```julia
m <= Setting(:saveroot, "path/to/my/save/root")
m <= Setting(:dataroot, "path/to/my/data/root")
```

Utility functions are provided to create paths to input/output files. These should be used
for best results.

```@docs
DSGE.inpath
DSGE.rawpath
DSGE.logpath
DSGE.workpath
DSGE.tablespath
DSGE.figurespath
```
