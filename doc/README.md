# DISCLAIMER

Copyright Federal Reserve Bank of New York.  You may reproduce, use, modify,
make derivative works of, and distribute and this code in whole or in part so
long as you keep this notice in the documentation associated with any
distributed works.   Neither the name of the Federal Reserve Bank of New York
(FRBNY) nor the names of any of the authors may be used to endorse or promote
works derived from this code without prior written permission.  Portions of the
code attributed to third parties are subject to applicable third party licenses
and rights.  By your use of this code you accept this license and any
applicable third party license.

OF ANY KIND, EITHER EXPRESS OR IMPLIED, INCLUDING WITHOUT LIMITATION ANY

FOR A PARTICULAR PURPOSE, EXCEPT TO THE EXTENT THAT THESE DISCLAIMERS ARE HELD
TO BE LEGALLY INVALID.  FRBNY IS NOT, UNDER ANY CIRCUMSTANCES, LIABLE TO YOU
FOR DAMAGES OF ANY KIND ARISING OUT OF OR IN CONNECTION WITH USE OF OR
INABILITY TO USE THE CODE, INCLUDING, BUT NOT LIMITED TO DIRECT, INDIRECT,
INCIDENTAL, CONSEQUENTIAL, PUNITIVE, SPECIAL OR EXEMPLARY DAMAGES, WHETHER

EQUITABLE THEORY, EVEN IF FRBNY HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH
DAMAGES OR LOSS AND REGARDLESS OF WHETHER SUCH DAMAGES OR LOSS IS FORESEEABLE.

# FRBNY DSGE Model (Version 990.2)

This module, written in Julia, estimates the model discussed in the
Liberty Street Economics blog post "The FRBNY DSGE Model Forecast."
This object-oriented, Julia-language implementation reproduces the
MATLAB code included in that post. Code to forecast the model is under
development and will be released upon completion.

# Model Design

The Julia implementation of the FRBNY model is designed around a
single object - the Model Object - which centralizes all information
about the model's parameters, states, equilibrium conditions, and
settings in a single data structure. The Model Object also keeps track
of file locations for all I/O operations. 

The general structure and basic functionality of the model object is
implemented by the AbstractDSGEModel type in
`abstractdsgemodel.jl`. For a specific model specification SPEC, the
model object is constructed by the function `ModelSPEC()` in
`mSPEC.jl`, which is a subtype of `AbstractDSGEModel`.

The following objects define a model:

- __Parameters__: Have values, bounds, fixed-or-not status, priors. An
  instance of the `Param` type houses all information about a given
  parameter in a single data structure.
- __States__: Mappings of names to indices (e.g. "π_t" -> 1)
- __Equilibrium Conditions__: A function that takes parameters and model
  indices, then returns Γ0, Γ1, Ψ, and Π (which fully describe the model in cannonical form)

These are enough to define the model structure. _Everything else_ is
essentially a function of these basics, and we can get to a forecast by
this chain:

- (Parameters + Model Indices + Eqcond Function) -> (TTT + RRR)
- (TTT + RRR + Data) -> Estimation
- (Estimation + TTT + RRR + Data) -> Forecast      (not yet implemented)


# Running the Code

## Running with Default Settings

So far, only the estimation step of the DSGE model has been
implemented. To run the estimation step from the Julia REPL, simply
create an instance of the model object and pass it to the `estimate`
function.

```julia
m = Model990()          # construct a model object
estimate(m)             # estimate the model
computeMoments(m)       # produce LaTeX tables of parameter moments
```

To change any of the model's default settings, parameters, equilibrium
conditions, etc., see "Implementation Details" for more specifics.


# Directory Structure

The directory structure follows Julia module conventions. In the
top-level directory (DSGE), you will find the folling subdirectory
tree (square brackets indicate future additions to the tree that will be added in future steps):

  - `doc/`: Documentation, including this README
  - `save/`: 
    - `input_data/`: This directory is referred to as `datapathroot` in the code.
      -`data/`:  Macroeconomic series formatted as an n x m Array{Float64,2}, where n is the number of observations and m is the number of series used as input.
      		 - `data_151030.h5`: input data vintage from October 30, 2015 (note that this is something of a misnomer because this data isnt from that date...)
	-`user/`: User-created files for model input. For instance, the user may specify a previously computed mode when `reoptimize(m)==false`, or a starting point for optimization when `reoptimize(m)==true`.
	  	   - `mode_in.h5`: Used as starting point for estimation when `reoptimize(m)==false`.
	     	   - `mode_in_optimized.h5`: Taken as the mode when `reoptimize(m)==true`.
	     	   - `hessian.h5`: Used as starting point for hessian calculation when `recalculate_hessian(m)==false`.
	     	   - `hessian_optimized.h5`: Taken as the hessian when `recalculate_hessian(m)==true`.
     - `output_data/`: 
         - `m990/`: Input/output files for the Model990 type. A model of type mSPEC will create its own save directory `mSPEC` at this  level in the directory tree.
              -`ss0/`: Subdirectory for subspec 0. We refer to this directory as `modelpathroot` in the code.
	      	  -`estimate`
			-`figures/`: Plots and other figures
			-`tables/`: LaTeX tables	          
			-`raw/`: Raw output from estimation step 
			     mode_out.h5: Optimized mode after running csminwel
			     hessian_out.h5: Hessian at the mode
			     sim_save.h5: Draws from posterior distribution			     
			-`work/`: HDF5 files created using `raw/` files as input
			     cov.h5: Covariance matrix for parameter draws from Metropolis-Hastings. Can be used as hessian matrix.
		  -[`forecast`]: Output for forecasts 
			-`figures/`: Plots and other figures
			-`tables/`: LaTeX tables	          
			-`raw/`: Raw output from forecast step
			-`work/`: HDF5 files created using `raw/` files as input
 		  -[`irfs`]: Impulse-response function outputs
		  -[`shockdecs`]: Shock decompositions
	      - [`ss1/`] Additional model subspecs will have subdirectories identical to `ss0` at this level in the directory tree. 
  - `src/`
     - `abstractdsgemodel.jl`: Defines the `AbstractDSGEModel` type.
     - `distributions_ext.jl`: Defines additional functions to return objects of type Distribution.
     - `DSGE.jl`: The main module file.
     - `estimate/`: Mode-finding and posterior sampling.
     - `models/`
           - `m990/`: Contains code to define and initialize version 990 of the FRBNY DSGE model.
              - `eqcond.jl`: Constructs `Model990`'s equilibrium condition matrices
              - `m990.jl`: Code for constructing a `Model990` object.
              - `measurement.jl`: Constructs `Model990`'s measurement equation matrices.
              - `subspecs.jl`: Code for model sub-specifications is defined here. See "Editing or Extending a Model" for details on constructing model sub-specifications.
           - [`m991/`]: Code for new subtypes of `AbstractDSGEModel` should be kept in directories at this level in the directory tree
     - `parameters.jl`: Implements the `AbstractParameter` type and its subtypes.
     - `settings.jl`: Implements the `Setting` type.
     - `solve/`: Solving the model; includes `gensys.jl` code.

  - `test/`: Module test suite.
   

# Implementation Details

This section describes important functions and implementation features in greater detail. If the user
is interested only in running the default model and reproducing the forecast
results, this section can be ignored.

This section focuses on what the code does and why, while the code itself
(including comments) provides detailed information regarding *how* these basic
procedures are implemented.

## The AbstractDSGEModel Type and the Model Object

The `AbstractDSGEModel` type provides a common infrastructure for all
model objects, which greatly facilitates the implementation of new
model specifications. Any concrete subtype of `AbstractDSGEModel` can
be passed to any function defined for `AbstractDSGEModel`, provided
that the concrete type has the fields that the function expects to be
available. 

`Model990` is one example of a concrete subtype of `AbstractDSGEModel`
that implements a single specification of the FRBNY DSGE
model. Infinitely more specifications are possible - see "Extending or
Editing a Model" below.

### Required Fields

#### Parameters and Steady-States
-`parameters::Vector{AbstractParameter}`: Vector of all time-invariant model parameters.
-`steady_state::Vector`: Model steady-state values, computed as a
function of elements of `parameters`.
-`keys::Dict{Symbol,Int}`: Maps human-readable names for all model
parameters and steady-states to their indices in `parameters` and
`steady_state`.

#### Inputs to the Measurement and Equilibrium Condition Equations
-`endogenous_states::Dict{Symbol,Int}`: Maps each state to a column
in the measurement and equilibrium condition matrices.
-`exogenous_shocks::Dict{Symbol,Int}`: Maps each shock to a column in
the measurement and equilibrium condition matrices.
-`expected_shocks::Dict{Symbol,Int}`: Maps each expected shock to a
column in the measurement and equilibrium condition matrices.
-`equilibrium_conditions::Dict{Symbol,Int}`: Maps each equlibrium
 condition to a row in the model's equilibrium condition matrices.
-`endogenous_states_postgensys::Dict{Symbol,Int}`: Maps lagged states
 to their columns in the measurement and equilibrium condition
 equations. These are added after Gensys solves the model.
-`observables::Dict{Symbol,Int}`: Maps each observable to a row in the
 model's measurement equation matrices.

#### Model Specification and Settings
-`spec::ASCIIString`: Model specification number (eg
 "m990"). Identifies a particular set of parameters, equilibrium
 conditions, and measurement equation (equivalently, a concrete model
 type - for example, models of type `Model990` should have spec = "m990".)
-`subspec::ASCIIString`: Model subspecification (eg "ss0"). Indicates any changes
 to parameter initialization from `spec`. See "Editing or Extending a
 Model" below for more details.
-`settings::Dict{Symbol,Setting}`: Settings/flags that affect
 computation without changing the economic or mathematical setup of
 the model.
-`test_settings::Dict{Symbol,Setting}`: Settings/flags for testing mode

#### Other Fields
-`rng::MersenneTwister`: Random number generator. By default, it is
 seeded to ensure replicability in algorithms that involve randomness
 (such as Metropolis-Hastings).
-`testing::Bool`: Indicates whether the model is in testing mode. If
 `true`, settings from `m.test_settings` are used in place of those in
 `m.settings`.
-`_filestrings::SortedDict{Symbol,AbstractString,ForwardOrdering}`:
An alphabetized list of setting identifier strings. These are
concatenated and appended to the filenames of all output files to
avoid overwriting the output of previous estimations/forecasts that
differ only in their settings, but not in their underlying
mathematical structure. See "Settings" below for more details.



## Defining Indices 

The model's equilibrium conditions and observables are represented as
fairly large matrices, and keeping track of which rows and columns
correspond to which states, shocks, equations, etc. can be
confusing. To improve clarity, we define several dictionaries
that map variable names to indices in these matrices:

- `endogenous_states`: Maps endogenous model states to columns in Γ0 and ZZ.
- `exogenous_shocks`:  Exogenous shocks
- `expected_shocks`:  Expectation shocks
- `equilibrium_conditions`: Equation indices
- `endogenous_states_postgensys`: Endogenous states, after model solution and
    system augmentation
- `observables`:  Indices of named observables to use in measurement equation

This approach has a number of advantages. Most importantly, it is
robust to inadvertent typos or indexing errors. Since the
actual index number doesn't matter to us, the user only needs to
define the names of their equilibrium conditions, states, and other
variables. Adding states is easy - we have only to add them to the
appropriate list in the model constructor, and they will be assigned
an index. No need to increment anything.

As an example, consider the model's equilibrium conditions. The
cannonical representation of the equilibrium conditions is 

`Γ0 s_t + Γ1 s_{t-1} + C + Ψ ε_t + Π η_t`

where `Γ0`, `Γ1`, `C`, `Ψ`, and `Π` are matrices of coefficients for `s_t`
(states at time $t$), `s_{t-1}` (lagged states), `ε_t` (iid shocks) and
`η_t` (expectational shocks). Each row of these matrices corresponds
to an equilibrium condition, which we define using a descriptive name
(for example, we name the consumption Euler equation `:euler`). States
(columns of Γ0 and Γ1), shocks (columns of Ψ), and expectational
states (columns Π) also have names. This is much more intuitive than
referring to the euler equation as `Γ0[1,:].`



## Parameters: The `AbstractParameter` Type

The `AbstractParameter` type implements our notion of a model
parameter: a time-invariant, unobserved value that has economic
significance in the model's equilibrium conditions. We estimate the
model to find the values of these parameters.

Though all parameters are time-invariant, each has different
features. Some parameters are scaled for use in the model's
equilibrium conditions and measurement equations.  During
optimization, parameters can be transformed from model space to the
real line via one of three different transformations. These
transformations are also defined as types, and require additional
information for each parameter. Finally, steady-state parameters are
not estimated directly, but are calculated as a function of other
parameters.

These various requirements are nicely addressed using a parameterized
type hierarchy. 

-`AbstractParameter{T<:Number}`: The common abstract supertype for all parameters.
    -`Parameter{T<:Number, U<:Transform}`: The abstract supertype for parameters that are directly estimated. 
        -`UnscaledParameter{T<:Number, U:<Transform}`: Concrete type for parameters that do not need to be scaled for equilibrium conditions.
        -`ScaledParameter{T<:Number, U:<Transform}`: Concrete type for parameters that are scaled for equilibrium conditions.
    -`SteadyStateParameter{T<:Number}`: Concrete type for steady-state parameters.


All `Parameter`s have the following fields:

- `key::Symbol`: Parameter name. For maximum clarity, `key` should
conform to the guidelines established in the DSGE Style Guide, in
CONTRIBUTING.md
- `value::T`: Parameter value. Initialized in model space (guaranteed
to be between `valuebounds`), but can be transformed between model
space and the real line via calls to `toreal` and `tomodel`.
- `valuebounds::Interval{T}`: Bounds for the parameter's value in
model space.
- `transform_parameterization::Interval{T}`: Parameters used to
transform `value` between model space and the real line.
- `transform::U`: Transformation used to transform `value` between
model space and real line.
- `prior::NullablePrior`: Prior distribution for parameter value.
- `fixed::Bool`: Indicates whether the parameter's value is fixed
rather than estimated.
- `description::AbstractString`: A short description of the
parameter's economic significance.
- `texLabel::AbstractString`: String for printing the parameter name
to LaTeX.

`ScaledParameters` also have the following fields:

-`scaledvalue::T`: Parameter value scaled for use in `eqcond.jl`
-`scaling::Function`: Function used to scale parameter value for use
in equilibrium conditions.

*Note:* Though not strictly necessary, defining a scaling with the
parameter object allows for much a much cleaner definition of the
equilibrium conditions.

Because the values of `SteadyStateParameter`s are directly computed as a
function of `ScaledParameter`s and `UnscaledParameter`s,
they only require 4 fields:

-`key::Symbol`
-`value::T`                    
-`description::AbstractString`
-`texLabel::AbstractString`



## `m.settings` and the `Setting` type

The `Setting` type implements computational settings that affect how
the code runs without affecting the mathematical definition of the
model. These include flags (e.g. whether or not to recompute the
hessian), parameterization for the Metropolis-Hastings algorithm
(e.g. number of times to draw from the posterior distribution), and
the vintage of data being used (`Setting` is a parametric type - a
`Setting{T<:Any}`, so Booleans, Numbers, and Strings can all be turned
into `Setting`s). They are stored centrally in the `settings`
dictionary within the model object.

Why implement a `Setting` type when we could put their values directly
into the source code or dictionary? The most obvious answer is that
the parametric type allows us to implement a single interface for all
`Setting`s (Booleans, Strings, etc), so that when we access a
particular setting during the estimation and forecast steps, we don't
have to think about the setting's type.

`Setting`s play an important role in addition to providing useful
abstraction. Estimating and forecasting the FRBNY DSGE model takes
many hours of computation time and creates a lot of output files. It
is useful to be able to compare model output from two different models
whose settings differ slightly (for example, consider two identical
models that use different vintages of data as input). A central
feature of the `Setting` type is a mechanism that generates unique,
meaningful filenames when code is executed with different
settings. Specifically, when a setting takes on a non-default value, a
user-defined setting code (along with the setting's value) are
appended to all output files generated during execution.

The `Setting{T<:Any}` type has the following fields:

- `key::Symbol`: Name of setting
- `value::T`: Value of setting 
- `savestring::Bool`: Indicates whether to append this setting's code
and value to output file names. If true, output file names will
include a suffix of the form _code1=val1_code2=val2_etc. where codes
are listed in alphabetical order.
- `code::AbstractString`: short string (<=4 characters) to print to output
file names when `savestring=true`.
- `description::AbstractString`: Short description of what the setting
is used for.

### Default Settings

#### I/O 

- `datapathroot::Setting{ASCIIString}`: The root directory for
model input data.
- `savepathroot::Setting{ASCIIString}`: The root directory for model output.
- `data_vintage`::Setting{ASCIIString}`: Data vintage identifier,
formatted YYMMDD (e.g. data from October 30, 2015 is identified by the
string "151030".) By default, `data_vintage` is set to the most recent
date of the files with name datapathroot/data/data_YYMMDD.h5. It is the only
setting printed to output filenames by default.

#### Anticipated Shocks
- `num_anticipated_shocks::Setting{Int}`: Number of anticipated policy shocks.
- `num_anticipated_shocks_padding::Setting{Int}`: Padding for
`num_anticipated_shocks`.
- `num_anticipated_lags::Setting{Int}`: Number of periods back to
incorporate zero bound expectations.
- `num_presample_periods::Setting{Int}`: Number of periods in the
presample

#### Estimation 
- `reoptimize::Setting{Bool}`: Whether to reoptimize the posterior
mode. If `false` (the default), `estimate()` reads in a previously
found mode.
- `recalculate_hessian::Setting{Bool}`: Whether to reecalculate the
hessian at the mode. If `false` (the default), `estimate()` reads in
a previously computed Hessian.

##### Metropolis-Hastings 
- `num_mh_simulations::Setting{Int}`: Number of draws from the
posterior distribution per block.
- `num_mh_blocks::Setting{Int}`: Number of blocks to run
Metropolis-Hastings.
- `num_mh_burn::Setting{Int}`: Number of blocks to discard as burn-in
for Metropolis-Hastings
- `mh_thinning_step::Setting{Int}`: Save every `mh_thinning_step`-th
draw in Metropolis-Hastings.


### Accessing Settings
The function `get_setting(m::AbstractDSGEModel, s::Symbol)` returns
the value of the setting `s` in `m.settings`. Some settings also
have explicit getter methods that take only the model object `m` as an argument:

*I/O settings*:
`modelpathroot(m)`
`datapathroot(m)`
`data_vintage(m)`

*Parallelization:*
`use_parallel_workers(m)`

*Estimation*
`reoptimize(m)`
`recalculate_hessian(m)`
`max_hessian_free_params(m)`

*Metropolis-Hastings*
`num_mh_blocks(m)`
`num_mh_simulations(m)`
`num_mh_burn(m)`
`mh_thinning_step(m)`	


### Overwriting Default Settings

To overwrite default settings added during model construction, a user
must define a new `Setting` object and overwrite the corresponding
entry in the model's `settings` dictionary using the `<=`
syntax. Individual fields of a pre-initialized setting object cannot
be modified. This immutability enforces the filenaming convention
described in the preceding paragraphs (the default parameters are
constructed without codes and are not printed to filename outputs to
avoid excessively long filenames). Therefore, we strongly suggest that
users who modify settings set `savestring=true` and define a
meaningful code when overwriting any default settings.

For example, overwriting `reoptimize` should look like this:
```julia
m = Model990()
# reoptimize(m) will return false by default
m <= Setting(:reoptimize, true, true, "reop", "whether to re-find the mode")
# reoptimize(m) will return true
```

### Test settings 

There are some settings that are always used when testing the
code. These are defined in `m.test_settings`, and returned by
`get_setting` when `m.testing=true` (if a setting is not in
`test_settings` but is referenced in the code, the value in `settings`
is used.)


## Estimation

**Main Function**: `estimate` in `src/estimate/estimate.jl`

**Purpose**: Finds modal parameter estimates and samples from posterior distribution. 

**Main Steps**: 

- *Initialization*: Read in and transform raw data from `save/input_data/`. 

- *Find Mode*: The main program will call the `csminwel` optimization
routine (located in `csminwel.jl`) to find modal parameter
estimates. Can optionally start estimation from a starting parameter
vector by specifying `save/input_data/mode_in.h5` If the
starting parameter vector is known to be optimized, the file should
be called `mode_in_optimized.h5`

- *Sample from Posterior*: Posterior sampling begins from the computed
mode, (or the provided mode if `reoptimize=false`), first computing
the Hessian matrix to scale the proposal distribution in the
Metropolis Hastings algorithm. Default settings for the number of
sampling blocks and the size of those blocks can be altered as
described in "Extending or Editing a Model".

**Remark**: In addition to saving each `mh_thinning_step`-th draw of the
parameter vector, the estimation program also saves the resulting
posterior value and transition equation matrices implied by each draw
of the parameter vector. This is to save time in the forecasting step
since that code can avoid recomputing those matrices. In addition, to
save space, all files in `save/input_data` and `save/output_data` are
HDF5 files.

# Extending or Editing a Model

Most users will want to extend or edit the DSGE model we provide in a
number of different ways. The most common changes we anticipate are
listed below, in decreasing order of complexity:

1. Add new parameters
2. Modify equilibrium conditions or measurement equations
3. Change the values of various parameter fields (i.e. initial `value`, `prior`, `transform`, etc)
4. Change the values of various computational settings (i.e. `reoptimize`, `num_mh_blocks`)

Points 1 and 2 often go together (adding a new parameter guarantees a
change in equilibrium conditions), and are such fundamental changes
that they increment the model specification nymber and require the
definition of a new subtype of `AbstractDSGEModel` (for instance,
`Model991`). See "Model specification" below for more details.

Any changes to the initialization of preexisting parameters are
defined as a new model *subspecification*, or *subspec*. While less
significant than a change to the model's equilibrium conditions,
changing the values of some parameter fields (especially priors) can
have economic significance over and above settings we use for
computational purposes. *For multiple reasons, parameter definitions
should not be modified in the model object's constructor.* First,
incrementing the model's subspecification number when parameters are
changed improves model-level (as opposed to code-level) version
control. Second, it avoids potential output filename collisions,
preventing the user from overwriting output from previous estimations
with the original parameters. The protocol for defining new
subspecifications is described below in "Model subspecifications".

Overriding default settings is described in the "Settings" section above.


## Model specification (`m.spec`) 

A particular model, which corresponds to a subtype of
`AbstractDSGEModel`, is defined as a set of parameters, equilibrium
conditions (defined by the `eqcond` function) and measurement
equations (defined by the `measurement` function). Therefore, the
addition of new parameters, states, or observables, or any changes to
the equilibrium conditions or measurement equations necessitate the
creation of a new subtype of `AbstractDSGEModel.`

To create a new model object, we recommend doing the following:

1. Duplicate the `m990` directory within the `src/models/`
directory. Name the new directory `mXXX.jl`, where `XXX` is your
chosen model specification number. Rename `m990.jl` in this directory to `mXXX.jl`.

2. In the `mXXX/` directory, change all references to `Model990` to `ModelXXX`.

3. Edit the `m990.jl`, `eqcond.jl`, and `measurement.jl` files as you
see fit. If adding new states, equilibrium conditions, shocks, or
observables, be sure to add them to the appropriate list in
`initialize_model_indices`.


## Model subspecifications (`m.subspec`)

Model990 subspecifications are initialized by overwriting initial
parameter definitions before the model object is fully
constructed. This happens via a call to `initialize_subspec` in the
`Model990` constructor. (Clearly, an identical protocol should be
followed for new model types as well.)

To create a new subspecification (e.g., subspec 1) of `Model990`, edit
the file `src/models/subspecs.jl` as follows (note that this example
is not actually subspecification 1 of Model990. In the source code,
our subspecification 5 is provided as additional example.):

1. Define a new function, `ss1`, that takes an object of type
`Model990` (not `AbstractDSGEModel`!) as an argument. In this
function, construct new parameter objects and overwrite existing model
parameters using the `<=` syntax. For example,

```julia
function ss1(m::Model990)

    m <= parameter(:ι_w, 0.000, (0.0, .9999), (0.0,0.9999), Untransformed(), Normal(0.0,1.0), fixed=false,
                   description="ι_w: Some description.",
                   texLabel="\\iota_w")

    m <= parameter(:ι_p, 0.0, fixed=true,
                   description= "ι_p: Some description"
                   texLabel="\\iota_p")

    
end
```


2. Add an `elseif` condition to `initialize_subspec`

```julia
    ...
    elseif subspec(m) == "ss1"   
        return ss1(m)               
    ...
```

To construct an instance of `Model990`, subspec `ss1`, call the
constructor for `Model990` with `ss1` as an
argument. For example,

```julia
m = Model990("ss1")
```

# Future Releases

This release marks only the beginning of our effort to convert all of
the MATLAB code released with the blog post referenced above. We plan
to release the following steps in the coming months:

1. Data step: Users can update the data, specify which vintage they
want to use in estimation/forecasting, and easily add new/different
observables to the dataset.

2. Forecast step: Given parameter draws, add code to compute forecasts,
shock decompositions, impulse responses, variance decompositions, and
counterfactuals. Add code to compute history of states and shocks.

3. Plotting step: Given forecasts other output from the Forecast step, generate the following figures:
   - Impulse-response functions
   - Shock decompositions
   - Shock histories
   - State histories

As these steps are under development, we would welcome improvements to
the existing code from the community. Some examples could be:
- Performance improvements
- Alternatives to algorithms used here (such as new optimization routines)


