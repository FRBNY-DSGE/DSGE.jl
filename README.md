# FRBNY DSGE Model (Version 990.2)
[![Build Status](https://travis-ci.org/FRBNY-DSGE/DSGE.jl.svg)](https://travis-ci.org/FRBNY-DSGE/DSGE.jl)

The *DSGE.jl* package implements the FRBNY DSGE model and provides general code
to estimate many user-specified DSGE models. The package is introduced in the
Liberty Street Economics blog post
[The FRBNY DSGE Model Meets Julia](http://libertystreeteconomics.newyorkfed.org/2015/12/the-frbny-dsge-model-meets-julia.html).

This Julia-language implementation mirrors the MATLAB code
included in the Liberty Street Economics blog post
[The FRBNY DSGE Model Forecast](http://libertystreeteconomics.newyorkfed.org/2015/05/the-frbny-dsge-model-forecast-april-2015.html).

FRBNY is currently working on extending the code to include forecasts and other
features. Extensions of the DSGE model code may be released in the future at
the discretion of FRBNY.

# Model Design

*DSGE.jl* is an object-oriented approach to solving the FRBNY DSGE model that
takes advantage of Julia's type system, multiple dispatch, package-handling
mechanism, and other features. A single model object centralizes all information
about the model's parameters, states, equilibrium conditions, and settings in a
single data structure. The model object also keeps track of file locations for
all I/O operations.

The following objects define a model:

- __Parameters__: Have values, bounds, fixed-or-not status, priors. An
  instance of the `AbstractParameter` type houses all information about a given
  parameter in a single data structure.
- __States__: Mappings of names to indices (e.g. `π_t` ➜ 1).
- __Equilibrium Conditions__: A function that takes parameters and model
  indices, then returns `Γ0`, `Γ1`, `C`, `Ψ`, and `Π` (which fully describe the
  model in canonical form).
- __Measurement Equation__: A function mapping states to observables.

These are enough to define the model structure. _Everything else_ is essentially
a function of these basics, and we can solve the model and forecast observables
via the following chain:

- Parameters + Model Indices + Equilibrium conditions ➜ Transition matrices
  in state-space form
- Transition matrices + Data ➜ Estimated parameter values
- Estimated parameters + Transition matrices + Data ➜ Forecast (not yet
  implemented)

# Running the Code

## Running with Default Settings

To run the estimation step in Julia, simply create an instance of the model object and pass
it to the `estimate` function --- see an [example](doc/examples/run_default.jl).

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

For more detail on changing the model's default settings, parameters, equilibrium
conditions, etc., see [Implementation Details](#implementation-details) for more specifics.

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

Note these locations can be overridden as desired:
```julia
m <= Setting(:saveroot, "path/to/my/save/root")
m <= Setting(:dataroot, "path/to/my/data/root")
```

# Input data

Given all of the hard work put into specifying the model, one should be able to maintain
the input data painlessly. To that extent, *DSGE.jl* provides facilities to download
appropriate vintages of data series from *FRED* (Federal Reserve Economic Data).

Note that a sample input dataset for use with model `m990` is provided; see [Sample input
data](#sample-input-data) for more details. To update this sample dataset for use with
model `m990`, see [Update sample input data](#update-sample-input-data).

## Setup

To take advantage of the ability to automatically download data series from FRED, set up
your FRED API access by following the directions
[here](https://github.com/micahjsmith/FredData.jl/blob/master/README.md).

## Loading data

At the most basic, loading data looks like this:

```
m = Model990()
df = load_data(m)
```

By default, `load_data` will look on the disk first to see if an appropriate vintage of data
is already present. If data on disk are not present, or if the data are invalid for any
reason, a fresh vintage will be downloaded from FRED and merged with the other data sources
specified. See `?load_data` for more details.

The resulting DataFrame `df` contains all the required data series for this model, fully
transformed. The first row is given by the Setting `date_presample_start` and the last row
is given by `date_mainsample_end`. The first `n_presample_periods` rows of `df` are the
presample.

Driver functions including `estimate` accept this `df` as an argument and convert it into a
`Matrix` suitable for computations using `df_to_matrix`, which sorts the data, ensures the
full sample is present, discards the date column, and sorts the observable columns according
to the `observables` field of the model object.

## Non-FRED data sources

Some data series may not be available from FRED or one may simply wish to use a different
data source, for whatever reason. The data sources and series are specified in the
`data_series` field of the model object. For each data source that is *not* `:fred`, a
well-formed CSV of the form `<source>_<yymmdd>.csv` is expected in the directory indicated
by `inpath(m, "data")`.  For example, the following might be the contents of a data source
for two series `:series1` and `:series2`:

```
date,series1,series2
1959-06-30,1.0,NaN
1959-09-30,1.1,0.5
etc.
```

Note that quarters are represented by the date of the *last* day of the quarter and missing
values are specified by `NaN`.

### Example

Let's consider an example dataset comprised of 10 macro series sourced from FRED and one
survey-based series sourced from, say, the Philadelphia Fed's [Survey of Professional
Forecasters](https://www.philadelphiafed.org/research-and-data/real-time-center/survey-of-professional-forecasters/historical-data/inflation-forecasts)
via Haver Analytics:

```
julia> m.data_series
Dict{Symbol,Array{Symbol,1}} with 2 entries:
 :spf   => [:ASACX10]
 :fred  => [:GDP, :PCE, ...] #etc
```

If the data vintage specified for the model is `151127` (Nov. 27, 2015), then the following
files are expected in `inpath(m, "data")`:

```
spf_151127.csv
fred_151127.csv
```

The FRED series will be downloaded and the `fred_151127.csv` file will be automatically
generated, but the `spf_151127.csv` file must be manually compiled as shown above:

```
date,ASACX10
1991-12-31,4.0
etc.
```

Now, suppose that we set the data vintage to `151222`, to incorporate the BEA's third
estimate of GDP. The `fred_151222.csv` file will be downloaded, but there are no updates to
the SPF dataset during this period. Regardless, the file `spf_151222.csv` must be present to
match the data vintage. The solution in this case is to manually copy and rename the older
SPF dataset. Although this is not an elegant approach, it is consistent with the concept of a
vintage as the data available at a certain point in time --- in this example, it just so
happens that the SPF data available on Nov. 27 and Dec. 22 are the same.

## Implementation

Let's quickly walk through the steps *DSGE* takes to create a suitable dataset.

First, a user provides a detailed specification of the data series and transformations used
for their model.
- the user specifies `m.observables`; the keys of this dictionary name the series to be used
    in estimating the model.
- the user specifies `m.data_series`; the keys of this dictionary name data sources, and the
    values of this dictionary are lists of mnemonics to be accessed from that data source.
    Note that these mnemonics do not correspond to observables one-to-one, but rather are
    usually series in *levels* that will be further transformed.
- the user specifies `m.data_transforms`; the keys of this dictionary name the series to be 
    constructed and match the keys of `m.observables` exactly; the values of this dictionary
    are functions that operate on a single argument (`levels`) which is a DataFrame of the
    series specified in `m.data_series`. These functions return a DataArray for a single
    series. These functions could do nothing (e.g. return `levels[:, :SERIES1]`) or
    perform a more complex transformation, such as converting to one quarter percent changes
    or adjusting into per-capita terms.
- the user adjusts data-related settings, such as `data_vintage`, `dataroot`,
    `date_presample_start`, `date_mainsample_end`, and `date_zlbregime_start`, and
    `use_population_forecast`.

Second, *DSGE* attempts to construct the dataset given this setup through a call to
`load_data`.
- *DSGE* checks the disk to see if a valid dataset is already stored. A dataset is valid if
    every series in `m.data_transforms` is present and the entire sample is contained (from
    `date_presample_start` to `date_mainsample_end`. If no valid dataset is already stored, the dataset will be recreated.
- In a preliminary stage, intermediate data series, as specified in `m.data_series`, are loaded in levels using
    `load_data_levels`. Note that these intermediate data series may contain more rows that
    the final dataset so that growth rates and other transformations can be successfully
    applied.
    - Data from the `:fred` data source are downloaded using `load_fred_data`. If an
        existing FRED vintage exists on disk, any required FRED series that is contained
        therein will be imported; any missing series will be downloaded directly from FRED
        using the *FredData* package.
    - Data from non-FRED data sources are read from disk, verified, and merged.
- The series in levels are transformed as specified in `m.data_transforms`.
    - To prepare for per-capita transformations, population data are filtered using
        `hpfilter`. Note that optionally, a population growth forecast from a service like
        Macroeconomic Advisers is appended to the recorded values before the filtering. Both
        filtered and unfiltered population levels and growth rates are added to the `levels`
        data frame.
    - The transformations are applied using the `levels` DataFrame as input.

Finally, the resulting dataset is saved to disk for future reference as `data_<yymmdd>.csv`
and the DataFrame is returned to the caller.
    
## Common pitfalls

Given the complexity of the data download, you may find that the dataset generated by
`load_data` is not exactly as you expect. Here are some common pitfalls to look out for:
- Ensure that the `data_vintage` model setting is as you expect. (Try checking
    `data_vintage(m)`.)
- If you are having a problem using *FredData*, ensure your API key is provided correctly
    and that there are no issues with your firewall, etc. Any issues with *FredData* proper
    should be reported on that project's page.
- Ensure that the `data_series` field of the model object is set as expected.
- Double check the transformations specified in the `data_transforms` field of the model
    object.
- Ensure that the keys of the `observables` and `data_transforms` fields of the model object
    match.
- Check the input files for [non-FRED data sources](#non-fred-data-sources). They should be
    in the directory indicated by `inpath(m, "data")`, be named appropriately given the
    vintage of data expected, and be formatted appropriately. One may have to copy and
    rename files of non-FRED data sources to match the specified vintage, even if the
    contents of the files would be identical.
- Look for any immediate issues in the final dataset saved (`data_<yymmdd>.csv`). If a data
    series in this file is all `NaN` values, then likely a non-FRED data source was not
    provided correctly.
- Ensure that the column names of the data CSV match the keys of the `observables` field of
    the model object.

## Sample input data

For more details on the sample input data provided -- which is used to estimate the provided
model `m990`, please see [Data](doc/Data.md).

For more details on using market interest rate expectations to treat the zero lower bound,
see [Anticipated Policy Shocks](doc/AnticipatedPolicyShocks.md). In particular, note that
our model, as used to compute the forecasts referenced in Liberty Street Economics posts,
is trained on data that includes six quarters of interest rate expectations. The user is
responsible for procuring interest rate expectations and appending it to the provided sample
data set, as discussed in the linked documentation here.

## Update sample input data

A sample dataset is provided for the 2015 Nov 27 vintage. To update this dataset:

1. See [above](#setup) to setup automatic data pulls using *FredData.jl*.
2. Specify the exact data vintage desired:
    ```
    julia> m <= Setting(:data_vintage, "yymmdd")
    ```
3. Create data files for the non-FRED data sources. For model `m990`, the required data
   files include `spf_<yymmdd>.csv` (with column `ASACX10`), `longrate_<yymmdd>.csv` (with
   column `FYCCZA`), and `fernald_<yymmdd>.csv` (with columns `TFPJQ` and `TFPKQ`). To
   include data on expected interest rates, the file `ois_<yymmdd>.csv` is also required. See
   [Data](doc/Data.md) for details on the series used and links to data sources.
4. Run `load_data(m)`; series from *FRED* will be downloaded and merged with the series from
   non-FRED data sources that you have already created. See [Common
   pitfalls](#common-pitfalls) for some potential issues.

# Implementation Details

This section describes important functions and implementation features in
greater detail. If the user is interested only in running the default model and
reproducing the forecast results, this section can be ignored.

This section focuses on what the code does and why, while the code itself
(including comments) provides detailed information regarding *how* these basic
procedures are implemented.

## Source Code Directory Structure

The source code directory structure follows Julia module conventions.

  - `doc/`: Code and model documentation
  - `src/`
     - `DSGE.jl`: The main module file.
     - `abstractdsgemodel.jl`: Defines the `AbstractModel` type.
     - `parameters.jl`: Implements the `AbstractParameter` type and its
       subtypes.
     - `settings.jl`: Implements the `Setting` type.
     - `distributions_ext.jl`: Defines additional functions to return objects of
       type Distribution.
     - `estimate/`: Mode-finding and posterior sampling.
     - [`xxx/`]: Other model functionality, such as forecasts, impulse response
       functions, and shock decompositions.
     - `models/`
           - `m990/`: Contains code to define and initialize version 990 of the
             FRBNY DSGE model.
              - `eqcond.jl`: Constructs `Model990` equilibrium condition
                matrices
              - `m990.jl`: Code for constructing a `Model990` object.
              - `measurement.jl`: Constructs `Model990` measurement equation
                matrices.
              - `subspecs.jl`: Code for model sub-specifications is defined
                here. See [Editing or Extending a Model](#editing-or-extending-a-model)
                for details on constructing model sub-specifications.
           - [`m991/`]: Code for new subtypes of `AbstractModel` should be kept
             in directories at this level in the directory tree
     - `solve/`: Solving the model; includes `gensys.jl` code.
  - `test/`: Module test suite.

## Reoptimizing

Generally, the user will want to reoptimize the parameter vector (and
consequently, calculate the Hessian at this new mode) every time they conduct
posterior sampling; that is, when:
- the input data are updated with a new quarter of observations or revised
- the model sub-specification is changed
- the model is derived from an existing model with different equilibrium
  conditions or measurement equation.

This behavior can be controlled more finely.

### Reoptimize from Starting Vector

Reoptimize the model starting from the parameter values supplied in use
in a specified file. Ensure that you supply an HDF5 file with a variable named
`params` that is the correct dimension and data type.
```julia
m = Model990()
params = load_parameters_from_file(m, "path/to/parameter/file.h5")
update!(m, params)
estimate(m)
```

### Skip Reoptimization Entirely

You can provide a modal parameter vector and optionally a Hessian matrix
calculated at that mode to skip the reoptimization entirely. These values are
usually computed by the user previously.

You can skip reoptimization of the parameter vector entirely.
```julia
m = Model990()
specify_mode!(m, "path/to/parameter/mode/file.h5")
estimate(m)
```
The `specify_mode!` function will update the parameter vector to the mode and
skip reoptimization. Ensure that you supply an HDF5 file with a variable named
`params` that is the correct dimension and data type. (See also the utility function
`load_parameters_from_file`.)

You can additionally skip calculation of the Hessian matrix entirely.
```julia
m = Model990()
specify_mode!(m, "path/to/parameter/mode/file.h5")
specify_hessian(m, "path/to/Hessian/matrix/file.h5")
estimate(m)
```
The `specify_hessian` function will cause `estimate` to read in the Hessian
matrix rather than calculating it directly.  Ensure that you supply an HDF5
file with a variable named `hessian` that is the correct dimension and data
type. Specifying the Hessian matrix but *not* the parameter mode results in
undefined behavior.

## The `AbstractModel` Type and the Model Object

The `AbstractModel` type provides a common interface for all model objects,
which greatly facilitates the implementation of new model specifications. Any
concrete subtype of `AbstractModel` can be passed to any function defined for
`AbstractModel`, provided that the concrete type has the fields that the
function expects to be available.

`Model990` is one example of a concrete subtype of `AbstractModel` that
implements a single specification of the FRBNY DSGE model.
See [Editing or Extending a Model](#editing-or-extending-a-model).

### Parameters and Steady-States
- `parameters::Vector{AbstractParameter}`: Vector of all time-invariant model
  parameters.
- `steady_state::Vector`: Model steady-state values, computed as a function of
  elements of `parameters`.
- `keys::Dict{Symbol,Int}`: Maps human-readable names for all model parameters
  and steady-states to their indices in `parameters` and
  `steady_state`.

### Inputs to the Measurement and Equilibrium Condition Equations
- `endogenous_states::Dict{Symbol,Int}`: Maps each state to a column in the
  measurement and equilibrium condition matrices.
- `exogenous_shocks::Dict{Symbol,Int}`: Maps each shock to a column in the
  measurement and equilibrium condition matrices.
- `expected_shocks::Dict{Symbol,Int}`: Maps each expected shock to a column in
  the measurement and equilibrium condition matrices.
- `equilibrium_conditions::Dict{Symbol,Int}`: Maps each equilibrium condition to
  a row in the model's equilibrium condition matrices.
- `endogenous_states_augmented::Dict{Symbol,Int}`: Maps lagged states to their
  columns in the measurement and equilibrium condition equations. These are
  added after `gensys` solves the model.
- `observables::Dict{Symbol,Int}`: Maps each observable to a row in the model's
  measurement equation matrices.

### Model Specification and Settings
- `spec::ASCIIString`: Model specification number (e.g. `"m990"`). Identifies a
  particular set of parameters, equilibrium conditions, and measurement equation
  (equivalently, a concrete model type - for example, models of type `Model990`
  would have `spec = "m990"`.)
- `subspec::ASCIIString`: Model sub-specification (e.g. `"ss0"`). Indicates any
  changes to parameter initialization from `spec`.
  See [Editing or Extending a Model](#editing-or-extending-a-model) for more
  details.
- `settings::Dict{Symbol,Setting}`: Settings/flags that affect computation
  without changing the economic or mathematical setup of the model.
- `test_settings::Dict{Symbol,Setting}`: Settings/flags for testing mode

### Other Fields
- `rng::MersenneTwister`: Random number generator. By default, it is
  seeded to ensure reproducibility in algorithms that involve randomness
  (such as Metropolis-Hastings).
- `testing::Bool`: Indicates whether the model is in testing mode. If `true`,
  settings from `m.test_settings` are used in place of those in `m.settings`.

## Defining Indices

The model's equilibrium conditions and observables are represented as fairly
large matrices, and keeping track of which rows and columns correspond to which
states, shocks, equations, etc. can be confusing. To improve clarity, we define
several dictionaries that map variable names to indices in these matrices:

- `endogenous_states`: Indices of endogenous model states
- `exogenous_shocks`: Indices of exogenous shocks
- `expected_shocks`: Indices of expectation shocks
- `equilibrium_conditions`: Indices of equilibrium condition equations
- `endogenous_states_augmented`: Indices of model states, after model solution
  and system augmentation
- `observables`:  Indices of named observables

This approach has a number of advantages. Most importantly, it is robust to
inadvertent typos or indexing errors. Since the actual index number doesn't
matter to us, the user only needs to define the names of their equilibrium
conditions, states, and other variables. Adding states is easy - we have only to
add them to the appropriate list in the model constructor, and they will be
assigned an index.

As an example, consider the model's equilibrium conditions. The canonical
representation of the equilibrium conditions is

```
Γ0 s_t = Γ1 s_{t-1} + C + Ψ ε_t + Π η_t
```

where `Γ0`, `Γ1`, `C`, `Ψ`, and `Π` are matrices of coefficients for `s_t`
(states at time `t`), `s_{t-1}` (lagged states), `ε_t` (exogenous shocks) and
`η_t` (expectational shocks). Each row of these matrices corresponds to an
equilibrium condition, which we define using a descriptive name (for example, we
name the consumption Euler equation `:euler`). States (columns of `Γ0` and
`Γ1`), exogenous shocks (columns of `Ψ`), and expectational shocks (columns
`Π`) also have names.

## Parameters: The `AbstractParameter` Type

The `AbstractParameter` type implements our notion of a model parameter: a
time-invariant, unobserved value that has economic significance in the model's
equilibrium conditions. We estimate the model to find the values of these
parameters.

Though all parameters are time-invariant, each has different features. Some
parameters are scaled for use in the model's equilibrium conditions and
measurement equations.  During optimization, parameters can be transformed from
model space to the real line via one of three different transformations. These
transformations are also defined as types, and require additional information
for each parameter. Finally, steady-state parameters are not estimated directly,
but are calculated as a function of other parameters.

These various requirements are nicely addressed using a parameterized type
hierarchy.

- `AbstractParameter{T<:Number}`: The common abstract supertype for all
  parameters.
    - `Parameter{T<:Number, U<:Transform}`: The abstract supertype for
      parameters that are directly estimated.
        - `UnscaledParameter{T<:Number, U:<Transform}`: Concrete type for
          parameters that do not need to be scaled for equilibrium conditions.
        - `ScaledParameter{T<:Number, U:<Transform}`: Concrete type for
          parameters that are scaled for equilibrium conditions.
    - `SteadyStateParameter{T<:Number}`: Concrete type for steady-state
      parameters.

All `Parameter`s have the following fields:

- `key::Symbol`: Parameter name. For maximum clarity, `key` should conform to
  the guidelines established in [CONTRIBUTING.md](CONTRIBUTING.md).
- `value::T`: Parameter value. Initialized in model space (guaranteed to be
  between `valuebounds`), but can be transformed between model space and the
  real line via calls to `transform_to_real_line` and
  `transform_to_model_space`.
- `valuebounds::Interval{T}`: Bounds for the parameter's value in model space.
- `transform_parameterization::Interval{T}`: Parameters used to transform
  `value` between model space and the real line.
- `transform::U`: Transformation used to transform `value` between model space
  and real line.
- `prior::NullablePrior`: Prior distribution for parameter value.
- `fixed::Bool`: Indicates whether the parameter's value is fixed rather than
  estimated.
- `description::AbstractString`: A short description of the parameter's economic
  significance.
- `tex_label::AbstractString`: String for printing the parameter name to LaTeX.

`ScaledParameters` also have the following fields:

- `scaledvalue::T`: Parameter value scaled for use in `eqcond.jl`
- `scaling::Function`: Function used to scale parameter value for use in
  equilibrium conditions.

*Note:* Though not strictly necessary, defining a scaling with the parameter
object allows for much a much cleaner definition of the equilibrium conditions.

Because the values of `SteadyStateParameter`s are directly computed as a
function of `ScaledParameter`s and `UnscaledParameter`s, they only require 4
fields:

- `key::Symbol`
- `value::T`
- `description::AbstractString`
- `tex_label::AbstractString`

## Model Settings

The `Setting` type implements computational settings that affect how the code
runs without affecting the mathematical definition of the model. These include
flags (e.g. whether or not to recompute the Hessian), parameterization for the
Metropolis-Hastings algorithm (e.g. number of times to draw from the posterior
distribution), and the vintage of data being used (`Setting` is a parametric
type - a `Setting{T<:Any}`, so Booleans, Numbers, and Strings can all be turned
into `Setting`s). They are stored centrally in the `settings` dictionary within
the model object.

Why implement a `Setting` type when we could put their values directly into the
source code or dictionary? The most obvious answer is that the parametric type
allows us to implement a single interface for all `Setting`s (Booleans, Strings,
etc.), so that when we access a particular setting during the estimation and
forecast steps, we don't have to think about the setting's type.

`Setting`s play an important role in addition to providing useful abstraction.
Estimating and forecasting the FRBNY DSGE model takes many hours of computation
time and creates a lot of output files. It is useful to be able to compare model
output from two different models whose settings differ slightly (for example,
consider two identical models that use different vintages of data as input). A
central feature of the `Setting` type is a mechanism that generates unique,
meaningful filenames when code is executed with different settings.
Specifically, when a setting takes on a non-default value, a user-defined
setting code (along with the setting's value) are appended to all output files
generated during execution.

The `Setting{T<:Any}` type has the following fields:

- `key::Symbol`: Name of setting
- `value::T`: Value of setting
- `print::Bool`: Indicates whether to append this setting's code and value to
  output file names. If true, output file names will include a suffix of the
  form `_key1=val1_key2=val2` etc. where codes are listed in alphabetical order.
- `code::AbstractString`: short string (4 characters or less) to print to output
  file names when `print=true`.
- `description::AbstractString`: Short description of what the setting is used
  for.

### Default Settings

See [defaults.jl](src/defaults.jl) for the complete description of default settings.

#### General

- `dataroot`: The root directory for
  model input data.
- `saveroot`: The root directory for model output.
- `use_parallel_workers`: Use available parallel workers in computaitons.
- `data_vintage`: Data vintage identifier, formatted
  `yymmdd`. By default, `data_vintage` is set to today's date. It is (currently) the only
  setting printed to output filenames by default.

#### Dates
- `date_presample_start`: Start date of pre-sample.
- `date_mainsample_start`: Start date of main sample.
- `date_zlbregime_start`: Start date of zero lower bound regime.
- `date_mainsample_end`: End date of main sample.
- `date_forecast_start`: Start date of forecast period.
- `date_forecast_end`: End date of forecast period.

#### Anticipated Shocks
- `n_anticipated_shocks`: Number of anticipated policy shocks.
- `n_anticipated_shocks_padding`: Padding for anticipated shocks.

#### Estimation
- `reoptimize`: Whether to reoptimize the posterior mode. If `true` (the default),
    `estimate()` begins reoptimizing from the model object's parameter vector.
- `calculate_hessian`: Whether to compute the Hessian. If `true` (the
    default), `estimate()` calculates the Hessian at the posterior mode.

#### Metropolis-Hastings
- `n_mh_simulations`: Number of draws from the posterior
    distribution per block.
- `n_mh_blocks`: Number of blocks to run Metropolis-Hastings.
- `n_mh_burn`: Number of blocks to discard as burn-in for
  Metropolis-Hastings.
- `mh_thin`: Metropolis-Hastings thinning step.

### Accessing Settings
The function `get_setting(m::AbstractModel, s::Symbol)` returns the value of the
setting `s` in `m.settings`. Some settings also have explicit getter methods
that take only the model object `m` as an argument. Note that not all are exported.

- I/O:
    - `saveroot(m)`,
    - `dataroot(m)`,
    - `data_vintage(m)`,
- Parallelization:
    - `use_parallel_workers(m)`
- Estimation:
    - `reoptimize(m)`,
    - `calculate_hessian(m)`,
- Metropolis-Hastings:
    - `n_mh_blocks(m)`,
    - `n_mh_simulations(m)`,
    - `n_mh_burn(m)`,
    - `mh_thin(m)`

### Overwriting Default Settings

To overwrite default settings added during model construction, a user must
define a new `Setting` object and update the corresponding entry in the
model's `settings` dictionary using the `<=` syntax. If the `print`, `code`, and
`description` fields of the new `Setting` object are not provided, the fields of the
existing setting will be maintained. If new values for `print`, `code`, and `description`
are specified, and if these new values are distinct from the defaults for those fields, the
fields of the existing setting will be updated.

For example, overwriting `use_parallel_workers` should look like this:
```julia
m = Model990()
m <= Setting(:use_parallel_workers, true)
```

## Estimation

Finds modal parameter values, calculate Hessian matrix at mode, and samples
from posterior distribution. See `estimate` in
[estimate.jl](src/estimate/estimate.jl).

**Main Steps**:

- *Initialization*: Read in and transform raw data from `save/input_data/`.

- *Reoptimize parameter vector*: The main program will call the `csminwel`
  optimization routine (located in `csminwel.jl`) to find modal parameter
  estimates.

- *Compute Hessian matrix*: Computing the Hessian matrix to scale the
  proposal distribution in the Metropolis-Hastings algorithm.

- *Sample from Posterior*: Posterior sampling is performed using the
  Metropolis-Hastings algorithm. A proposal distribution is constructed centered
  at the posterior mode and with proposal covariance scaled by the inverse of
  the Hessian matrix. Settings for the number of sampling blocks and the size of
  those blocks can be altered as described in
  [Editing or Extending a Model](#editing-or-extending-a-model).

**Remark**: In addition to saving each `mh_thin`-th draw of the parameter
vector, the estimation program also saves the resulting posterior value and
transition equation matrices implied by each draw of the parameter vector. This
is to save time in the forecasting step since that code can avoid recomputing
those matrices.

# Editing or Extending a Model

Users may want to extend or edit `Model990` in a number of different ways.
The most common changes are listed below, in decreasing order of complexity:

1. Add new parameters
2. Modify equilibrium conditions or measurement equations
3. Change the values of various parameter fields (i.e. initial `value`,
   `prior`, `transform`, etc.)
4. Change the values of various computational settings (i.e. `reoptimize`,
   `n_mh_blocks`)

Points 1 and 2 often go together (adding a new parameter guarantees a change in
equilibrium conditions), and are such fundamental changes that they increment
the model specification number and require the definition of a new subtype of
`AbstractModel` (for instance, `Model991`).
See [Model specification](#model-specification-mspec) for more details.

Any changes to the initialization of preexisting parameters are defined as a new
model *sub-specification*, or *subspec*. While less significant than a change to
the model's equilibrium conditions, changing the values of some parameter fields
(especially priors) can have economic significance over and above settings we
use for computational purposes. **Parameter definitions should not be modified
in the model object's constructor.** First, incrementing the model's
sub-specification number when parameters are changed improves model-level (as
opposed to code-level) version control. Second, it avoids potential output
filename collisions, preventing the user from overwriting output from previous
estimations with the original parameters. The protocol for defining new
sub-specifications is described in
[Model sub-specifications](#model-sub-specifications-msubspec).

Overriding default settings is described in the [Settings](#model-settings)
section above.

## Model specification (`m.spec`)

A particular model, which corresponds to a subtype of `AbstractModel`, is
defined as a set of parameters, equilibrium conditions (defined by the `eqcond`
function) and measurement equations (defined by the `measurement` function).
Therefore, the addition of new parameters, states, or observables, or any
changes to the equilibrium conditions or measurement equations necessitate the
creation of a new subtype of `AbstractModel.`

To create a new model object, we recommend doing the following:

1. Duplicate the `m990` directory within the [models](src/models/) directory. Name the
new directory `mXXX.jl`, where `XXX` is your chosen model specification number
or string.  Rename `m990.jl` in this directory to `mXXX.jl`.

2. In the `mXXX/` directory, change all references to `Model990` to `ModelXXX`.

3. Edit the `m990.jl`, `eqcond.jl`, and `measurement.jl` files as you see fit.
If adding new states, equilibrium conditions, shocks, or observables, be sure to
add them to the appropriate list in `init_model_indices`.

4. Open the module file, `src/DSGE.jl`. Add `ModelXXX` to the list of functions
to export, and include each of the files in `src/model/mXXX`.

## Model sub-specifications (`m.subspec`)

`Model990` sub-specifications are initialized by overwriting initial parameter
definitions before the model object is fully constructed. This happens via a
call to `init_subspec` in the `Model990` constructor. (Clearly, an identical
protocol should be followed for new model types as well.)

To create a new sub-specification (e.g., subspec 1) of `Model990`, edit the file
`src/models/subspecs.jl` as follows (note that this example is not actually
sub-specification `1` of `Model990`. In the source code, our sub-specification
`5` is provided as additional example.):

1. Define a new function, `ss1`, that takes an object of type `Model990` (not
`AbstractModel`!) as an argument. In this function, construct new parameter
objects and overwrite existing model parameters using the `<=` syntax. For
example,

```julia
function ss1(m::Model990)
    m <= parameter(:ι_w, 0.000, (0.0, .9999), (0.0,0.9999), DSGE.Untransformed(), Normal(0.0,1.0), fixed=false,
                   description="ι_w: Some description.",
                   tex_label="\\iota_w")
    m <= parameter(:ι_p, 0.0, fixed=true,
                   description= "ι_p: Some description"
                   tex_label="\\iota_p")
end
```

2. Add an `elseif` condition to `init_subspec`:

```julia
    ...
    elseif subspec(m) == "ss1"
        return ss1(m)
    ...
```

To construct an instance of `Model990`, `ss1`, call the constructor
for `Model990` with `ss1` as an argument. For example,

```julia
m = Model990("ss1")
```

# Acknowledgements
Developers of this package at [FRBNY](https://www.newyorkfed.org/research)
include

* [Pearl Li](https://github.com/pearlzli)
* [Erica Moszkowski](https://github.com/emoszkowski)
* [Micah Smith](https://github.com/micahjsmith)

Contributors to this package at [QuantEcon](http://quantecon.org) include

* [Zac Cranko](https://github.com/ZacCranko)
* [Spencer Lyon](https://github.com/spencerlyon2)
* [Pablo Winant](http://www.mosphere.fr/)

The `gensys` and `csminwel` routines in [gensys.jl](src/solve/gensys.jl) and
[csminwel.jl](src/estimate/csminwel.jl) are based on routines originally
copyright [Chris Sims](http://www.princeton.edu/~sims). The files are released
here with permission of Chris Sims under the BSD-3 [license](LICENSE).

The `kalman_filter` routine is loosely based on a version of the
Kalman filter algorithm originally copyright Federal Reserve Bank of Atlanta
and written by [Iskander Karibzhanov](http://karibzhanov.com). The files are
released here with permission of the Federal Reserve Bank of Atlanta under the
BSD-3 [license](LICENSE).
