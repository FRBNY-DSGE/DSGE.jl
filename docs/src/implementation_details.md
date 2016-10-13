# Implementation Details

This section describes important functions and implementation features in
greater detail. If the user is interested only in running the default model and
reproducing the estimation results, this section can be ignored. Additional documentation
can also be found in function documentation or in-line.

This section focuses on what the code does and why, while the code itself
(including comments) provides detailed information regarding *how* these basic
procedures are implemented.

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

