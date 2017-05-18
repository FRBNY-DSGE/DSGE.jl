# Implementation Details

```@meta
CurrentModule = DSGE
```

This section describes important functions and implementation features in
greater detail. If the user is interested only in running the default model and
reproducing the estimation results, this section can be ignored. Additional documentation
can also be found in function documentation or in-line.

This section focuses on what the code does and why. Docstrings and the code itself
(including comments) provide detailed information regarding *how* these basic
procedures are implemented.

## The `AbstractModel` Type

The `AbstractModel` type provides a common interface for all model objects,
which greatly facilitates the implementation of new model specifications. Any
concrete subtype of `AbstractModel` can be passed to any function defined for
`AbstractModel`, provided that the concrete type has the fields that the
function expects to be available.

`Model990` is one example of a concrete subtype of `AbstractModel` that
implements a single specification of the New York Fed DSGE model. All model
objects must have these fields so that the interface for `AbstractModel` objects
works correctly.  See [Editing or Extending a Model](@ref
editing-extending-model) for more detail.

```@docs
Model990
```

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

## The `AbstractParameter` Type

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

All `Parameter`s have the fields defined in `UnscaledParameter`:

```@docs
UnscaledParameter
```

`ScaledParameters` also have the following fields:

- `scaledvalue::T`: Parameter value scaled for use in `eqcond.jl`
- `scaling::Function`: Function used to scale parameter value for use in
  equilibrium conditions.

*Note:* Though not strictly necessary, defining a scaling with the parameter
object allows for much a much cleaner definition of the equilibrium conditions.

Because the values of `SteadyStateParameter`s are directly computed as a
function of `ScaledParameter`s and `UnscaledParameter`s, they only require 4
fields:

```@docs
SteadyStateParameter
```


## The `Observable` and `PseudoObservable` Types

We similarly encapsulate information about observables and pseudo-observables
(unobserved linear combinations of states, e.g. the output gap) into the
`Observable` and `PseudoObservable` types. Each type has identifier fields
`key`, `name`, and `longname`.

Most importantly, both `Observable`s and `PseudoObservable`s include the
information needed for transformations to and from model units. For
`Observable`s, these are the `input_series`, `fwd_transform`, and
`rev_transform` fields. "Forward transformations" are applied to transform
the raw input data series specified in `input_series` to model units. The
model is estimated and forecasted in model units, and then we apply "reverse
transformations" to get human-readable units before computing means and bands or
plotting. Pseudo-observables are not observed, so they do not have
`input_series` or `fwd_transform`s, but they may however have `rev_transform`s.

As an example, the `:obs_gdp` `Observable` uses as `input_series` aggregate
nominal GDP in levels, the GDP price index, and population in levels, all from
FRED. These series are `fwd_transform`ed to get quarter-over-quarter log growth
rates of per-capita real GDP, which are the `Observable`'s model units. The
reverse transformation then converts `:obs_gdp` into annualized
quarter-over-quarter percent changes of *aggregate* real GDP.

```@docs
Observable
PseudoObservable
```


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
Estimating and forecasting the New York Fed DSGE model takes many hours of
computation time and creates a lot of output files. It is useful to be able to
compare model output from two different models whose settings differ slightly
(for example, consider two identical models that use different vintages of data
as input). A central feature of the `Setting` type is a mechanism that generates
unique, meaningful filenames when code is executed with different settings.
Specifically, when a setting takes on a non-default value, a user-defined
setting code (along with the setting's value) are appended to all output files
generated during execution.

The `Setting{T<:Any}` type is defined as follows:

```@docs
Setting
```

To update the value of an existing function, the user has two
options. First, the user may use the `<=` syntax as shown in
the [Running with Default Settings](@ref) section. However, for this
to work properly, it is essential that the setting's `key` field be
exactly the same as that of an existing entry in
`m.settings`. Otherwise, an additional entry will be added to
`m.settings` and the old setting will be the one accessed from other
all routines. A potentially safer, though clunkier, option is to use the [`update!`](@ref) method.

## Type Interfaces

### `AbstractModel` Interface

```@docs
DSGE.update!
DSGE.transform_to_model_space!
DSGE.load_parameters_from_file
DSGE.specify_mode!
DSGE.specify_hessian
```

### `Parameter` Interface

```@autodocs
Modules = [DSGE]
Pages = ["parameters.jl"]
Order = [:function]
```

### `Setting` Interface

```@autodocs
Modules = [DSGE]
Pages = ["settings.jl"]
Order = [:function]
```
