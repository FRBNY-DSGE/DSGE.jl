# [Estimation](@id estimation-step)

```@meta
CurrentModule = DSGE
```

## Procedure

The goal of the estimation step is to sample from the posterior
distribution of the model parameters. DSGE.jl provides two estimation routines.
The default is a Metropolis-Hastings
sampler to do this, which requires as a proposal covariance matrix and the
Hessian matrix corresponding to the posterior mode.
The second routine is a Sequential Monte Carlo sampler, which is called
from the [SMC.jl](https://github.com/FRBNY-DSGE/SMC.jl) package.
Both routines implement adaptive proposal densities and parameter blocking.[^1]

The function `estimate` implements the entire procedure for either routine.
Below, we explain the MH algorithm.
For documentation of the SMC algorithm, see [here](https://frbny-dsge.github.io/SMC.jl/latest/).

**Main Steps**:

- *Initialization*: Read in and transform raw data from `save/input_data/`. See [Input Data](@ref input-data-step) for more details.

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

**Remark**: For the MH sampler, in addition to saving each `mh_thin`-th draw of the parameter
vector, the estimation program also saves the resulting posterior value and
transition equation matrices implied by each draw of the parameter vector. This
is to save time in the forecasting step since that code can avoid recomputing
those matrices.

For the SMC sampler, we save a `jld2` file containing a `Cloud` object which holds
any relevant information about the particles approximating the posterior,
the normalized weights of the particles, and the unnormalized weights of particles.
We also save a draw of parameters to an `h5` file.


To run the entire procedure, the user simply calls the `estimate` routine:

```@docs
DSGE.estimate
```

## Metropolis-Hastings Sampler

### Computing the Posterior

In DSGE.jl, the function `posterior` computes the value of the posterior
distribution at a given parameter vector. It calls the `likelihood` function,
which in turn calls the `filter` routine. See [Estimation routines](@ref) for
more details on these functions.

We implement the Kalman Filter via the `filter` function to compute the
log-likelihood, and add this to the log prior to obtain the log posterior. See
[StateSpaceRoutines.jl](https://github.com/FRBNY-DSGE/StateSpaceRoutines.jl) for
a model-independent implementation of the Kalman filter.


### [Optimizing or Reoptimizing](@id estimation-reoptimizing)

Generally, the user will want to reoptimize the parameter vector (and consequently,
calculate the Hessian at this new mode) every time they conduct posterior sampling; that is,
when:

- the input data are updated with a new quarter of observations or revised
- the model sub-specification is changed
- the model is derived from an existing model with different equilibrium conditions or
  measurement equation.

This behavior can be controlled more finely.

#### Reoptimize from a Specified Starting Vector

Reoptimize the model starting from the parameter values supplied in a specified file.
Ensure that you supply an HDF5 file with a variable named `params` that is the correct
dimension and data type.

```julia
m = Model990()
params = load_parameters_from_file(m, "path/to/parameter/file.h5")
update!(m, params)
estimate(m)
```

#### Skip Reoptimization Entirely

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

#### Random Walk Metropolis-Hastings

For relatively simple problems, a random walk MH sampler is sufficient and avoids
unnecessary computations like calculating the Hessian. Suppose we want to estimate
all the parameters of a DSGE model `m`. The following code
implements RWMH for `m`.

```julia
m = Model990()
m <= Setting(:hessian_path, "path/to//matrix/with/right/dimensions/saved/as/mh_hessian.h5")
estimate(m; proposal_covariance = Matrix{Float64}(I,size(m.parameters)))
```

The saved Hessian needs to have the same dimensions as the number of parameters. The simplest
option is save an identity matrix. The proposal covariance also needs to have the
same dimensions as the number of parameters. If the user does not want to
estimate every parameter (i.e. some parameters are fixed), then the user
needs to zero out the rows of the proposal covariance that correspond to fixed parameters.


### Calculating the Hessian

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

### Estimation routines

#### Prior, Likelihood and Posterior calculations

```@docs
DSGE.prior
DSGE.likelihood
DSGE.posterior
DSGE.posterior!
```

#### Optimization

See [Optimization](@ref algs-optimization)

#### Full Estimation Routine

See [`estimate`](@ref)

#### Output Analysis

```@autodocs
Modules = [DSGE]
Pages   = ["moments.jl"]
Order   = [:function, :type]
```

## SMC Sampler

See [here](https://frbny-dsge.github.io/SMC.jl/latest/) for the settings
adjusting the SMC algorithm. To use these with a DSGE model object,
either add them after the definition of a model object
as a `Setting` or add them directly in the definition of the model type. For example,
the following code sets the sampler as SMC
and the number of particles used by SMC as 10,000:

```julia
m = Model990()
m <= Setting(:sampling_method, :SMC)
m <= Setting(:n_particles, 10000)
```

[^1]: We document the details of implementing adaptive proposal densities and
      parameter blocking in [SMC.jl](https://github.com/FRBNY-DSGE/SMC.jl).