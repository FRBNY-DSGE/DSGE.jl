# [Estimation](@id estimation-step)

```@meta
CurrentModule = DSGE
```

## Procedure

The goal of the estimation step is to sample from the posterior
distribution of the model parameters. DSGE.jl uses either a Metropolis-Hastings or Sequential Monte Carlo
sampler to do this. The Metropolis-Hastings sampler requires as a  proposal covariance matrix the
Hessian matrix corresponding to the posterior mode. Depending on the value of the model setting `sampling_method` (which can be checked via `get_setting(m, :sampling_method)` and set via `m <= Setting(:sampling_method, :SMC)`, the function `estimate` will either find the mode, compute the Hessian, and run Metropolis-Hastings (if `get_setting(m, :sampling_method)==:MH)`) or run Sequential Monte Carlo (if `get-setting(m, :sampling_method)==:SMC`).

### Metropolis Hastings

**Main Steps**:

- *Initialization*: Read in and transform raw data from `save/input_data/`. See [Input Data](@ref) for more details.

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
  [Editing or Extending a Model](@ref editing-extending-model).

**Remark**: In addition to saving each `mh_thin`-th draw of the parameter
vector, the estimation program also saves the resulting posterior value and
transition equation matrices implied by each draw of the parameter vector. This
is to save time in the forecasting step since that code can avoid recomputing
those matrices.

To run the entire procedure, the user simply calls the `estimate` routine:

```@docs
DSGE.estimate
```

#### [Optimizing or Reoptimizing](@id estimation-reoptimizing)

Generally, the user will want to reoptimize the parameter vector (and consequently,
calculate the Hessian at this new mode) every time they conduct posterior sampling; that is,
when:

- the input data are updated with a new quarter of observations or revised
- the model sub-specification is changed
- the model is derived from an existing model with different equilibrium conditions or
  measurement equation.

This behavior can be controlled more finely.

##### Reoptimize from a Specified Starting Vector

Reoptimize the model starting from the parameter values supplied in a specified file.
Ensure that you supply an HDF5 file with a variable named `params` that is the correct
dimension and data type.

```julia
m = Model990()
params = load_parameters_from_file(m, "path/to/parameter/file.h5")
update!(m, params)
estimate(m)
```

##### Skip Reoptimization Entirely

You can provide a modal parameter vector and optionally a Hessian matrix calculated at that
mode to skip the reoptimization entirely. These values are usually computed by the user
previously.

You can skip reoptimization of the parameter vector entirely.

```julia
m = Model990()
specify_mode!(m, "path/to/parameter/mode/file.h5")
estimate(m)
```

The `specify_mode!` function will update the parameter vector to the mode and skip
reoptimization by setting the `reoptimize` model setting. Ensure that you supply an HDF5
file with a variable named `params` that is the correct dimension and data type. (See also
the utility function `load_parameters_from_file`.)

#### Calculating the Hessian

By default, `estimate` will recompute the Hessian matrix. You can skip
calculation of the Hessian matrix entirely if you provide a file with
a Hessian that has been pre-computed.

```julia
m = Model990()
specify_mode!(m, "path/to/parameter/mode/file.h5")
specify_hessian(m, "path/to/Hessian/matrix/file.h5")
estimate(m)
```

The `specify_hessian` function will cause `estimate` to read in the Hessian matrix rather
than calculating it directly.  Ensure that you supply an HDF5 file with a variable named
`hessian` that is the correct dimension and data type. Specifying the Hessian matrix but
*not* the parameter mode results in undefined behavior.

See [Hessian Approximation] for more details on the Hessian computation.


#### Computing the Posterior

In DSGE.jl, the function `posterior` computes the value of the posterior
distribution at a given parameter vector. It calls the `likelihood` function,
which in turn calls the `filter` routine. See [Estimation routines](@ref) for
more details on these functions.

We implement the Kalman Filter via the `filter` function to compute the
log-likelihood, and add this to the log prior to obtain the log posterior. See
[StateSpaceRoutines.jl](https://github.com/FRBNY-DSGE/StateSpaceRoutines.jl) for
a model-independent implementation of the Kalman filter.

### Sequential Monte Carlo

#### Procedure

**Main Steps**:

- *Initialization*: Read in and transform raw data from `save/input_data/`. See [Input Data](@ref) for more details.

- *Sample from Posterior*: Posterior sampling is performed using the Sequential Monte Carlo algorithm. A cloud of particles and their weights is initialized from the prior. These particles are mutated and their weights adjusted based on each value's likelihood. Settings for the number of particles, tempering schedule, and various other settings can be altered as described in [Editing or Extending a Model](@ref editing-extending-model).

To run the entire procedure, the user simply calls the `estimate` routine:

```@docs
DSGE.estimate
```

The user can also call `smc(m, df)` explicitly (where `m` is a model object and `df` is a `DataFrame`

#### Computing the Posterior

In DSGE.jl, the function `posterior` computes the value of the posterior
distribution at a given parameter vector. It calls the `likelihood` function,
which in turn calls the `filter` routine. See [Estimation routines](@ref) for
more details on these functions.

We implement the Kalman Filter via the `filter` function to compute the
log-likelihood, and add this to the log prior to obtain the log posterior. See
[StateSpaceRoutines.jl](https://github.com/FRBNY-DSGE/StateSpaceRoutines.jl) for
a model-independent implementation of the Kalman filter.


## Estimation routines

### Prior, Likelihood and Posterior calculations

```@docs
DSGE.prior
DSGE.likelihood
DSGE.posterior
DSGE.posterior!
```

### Optimization

See [Optimization](@ref algs-optimization)

### Full Estimation Routine

See [`estimate`](@ref)

### Output Analysis

```@autodocs
Modules = [DSGE]
Pages   = ["moments.jl"]
Order   = [:function, :type]
```
