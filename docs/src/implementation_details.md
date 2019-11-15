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

As of DSGE.jl v0.7.3, many types housed in the DSGE.jl package have been moved to
[ModelConstructors.jl](https://github.com/FRBNY-DSGE/ModelConstructors.jl).
The following types now belong in ModelConstructors.jl:

- `AbstractModel`
- `AbstractParameter`
    - `Parameter{T<:Number, U<:Transform}`: The abstract supertype for
      parameters that are directly estimated.
        - `UnscaledParameter{T<:Number, U:<Transform}`: Concrete type for
          parameters that do not need to be scaled for equilibrium conditions.
        - `ScaledParameter{T<:Number, U:<Transform}`: Concrete type for
          parameters that are scaled for equilibrium conditions.
    - `SteadyStateParameter{T<:Number}`: Concrete type for steady-state
      parameters.
- `Setting`
- `Observable`
- `PseudoObservable`
- Types and functions used to define and work with priors

We refer users to the documentation provided for ModelConstructors.jl
for information about the implementation of these types. Below, we
document the implementation of DSGE.jl specific types.


## The `AbstractDSGEModel` Type

The `AbstractModel` type has been rewritten to be an abstract type
for any model with parameters. We have replaced
the `AbstractModel` type in DSGE.jl with the `AbstractDSGEModel` type, which
is a subtype of `AbstractModel` and includes various methods that standard
DSGE models need.

The `AbstractDSGEModel` type provides a common interface for all DSGE model objects,
which greatly facilitates the implementation of new model specifications. Any
concrete subtype of `AbstractDSGEModel` can be passed to any function defined for
`AbstractDSGEModel`, provided that the concrete type has the fields that the
function expects to be available.

`Model990` is one example of a concrete subtype of `AbstractDSGEModel` that
implements a single specification of the New York Fed DSGE model. All model
objects must have these fields so that the interface for `AbstractDSGEModel` objects
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


## Type Interfaces

### `AbstractDSGEModel` Interface

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
