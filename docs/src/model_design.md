# Model Design

*DSGE.jl* is an object-oriented approach to solving the New York Fed DSGE model
that takes advantage of Julia's type system, multiple dispatch, package-handling
mechanism, and other features. A single model object centralizes all information
about the model's parameters, states, equilibrium conditions, and settings in a
single data structure. The model object also keeps track of file locations for
all I/O operations.

The following objects define a model:

- **Parameters**: Have values, bounds, fixed-or-not status, priors. An
  instance of the `AbstractParameter` type houses all information about a given
  parameter in a single data structure. See
  [The `AbstractParameter` Type](@ref).
- **Model Indices**: Mappings of state, shock, observable, and pseudo-observable
  names to indices (e.g. `y_t` -> 1). See [Defining Indices](@ref).
- **Observables and PseudoObservables**: Mapping of names to indices, as well as
  information necessary for transformations. See
  [The `Observable` and `PseudoObservable` Types](@ref).
- **Equilibrium Conditions**: A function that takes parameters and model
  indices, then returns `Γ0`, `Γ1`, `C`, `Ψ`, and `Π` (which fully describe the
  model in canonical form).
- **Measurement Equation**: A function mapping states to observables.
- **Pseudo-Measurement Equation**: A function mapping states to what we call
  "pseudo-observables", i.e. linear combinations of existing states which are
  not observed. (Note that this is not strictly required to implement a model,
  but we often use the pseudo-measurement equation instead of adding new states
  in order to achieve a more parsimonious model.)

These are enough to define the model structure. *Everything else* is essentially
a function of these basics, and we can solve the model and forecast observables
via the following chain:

- Parameters + Model Indices + Equilibrium conditions -> Transition matrices
  in state-space form
- Transition matrices + Data -> Estimated parameter values
- Estimated parameters + Transition matrices + Data -> Forecast
