# Special Model Types

```@meta
CurrentModule = DSGE
```

In addition to the `AbstractDSGEModel` type, *DSGE.jl* has several special types that implement
models which are designed to interface with DSGEs.

## The `PoolModel` Type
Unlike the other models contained in DSGE, the `PoolModel` type is not a proper DSGE model.
It is a wrapper object for different methods to perform model averaging for two different models,
which do not have to be DSGE models.
For example, a user could average two different vector auto-regressions.
Generally, a user only needs to provide the predictive density scores
of the two models that the user wants to average. The reason is that we
treat the predictive density scores as non-FRED observables. This approach
makes interfacing with the rest of the machinery provided by *DSGE.jl* very simple.

```@docs
PoolModel
```

See [Del Negro et al. (2016)](https://www.sciencedirect.com/science/article/pii/S0304407616300094#f000005) for theoretical details on the model averaging methods listed in the documentation.

To facilitate analysis with the `PoolModel` type, we also provide the following functions.
```@docs
estimate_bma
sample_λ
propagate_λ
compute_Eλ
```

## DSGE-VARs and the `DSGEVAR` Type

We can approximate the dynamics of a linearized DSGE with a VAR(``p``), where ``p`` is the
number of lags. This approximation gives a mapping from the parameters of a DSGE
to the parameters of a VAR (the coefficients and innovations variance-covariance matrix).
Since the number of parameters in a VAR are generally larger than the number of
parameters in a DSGE, this mapping can be interpreted as cross-restrictions
imposed by a DSGE on the parameters of a VAR. A DSGE-VAR combines a DSGE with a VAR
to, among other reasons, evaluate the mis-specification of the DSGE and improve
the DSGE's forecasting performance.

For more details on the theory, see
[Del Negro and Schorfheide (2004)](https://onlinelibrary.wiley.com/doi/full/10.1111/j.1468-2354.2004.00139.x),
[Del Negro and Schorfheide (2006)](https://www.newyorkfed.org/medialibrary/media/research/economists/delnegro/erq206_delnegro.pdf), and
[Del Negro and Schorfheide (2009)](https://www.aeaweb.org/articles?id=10.1257/aer.99.4.1415).

We implement DSGE-VARs with the `DSGEVAR` concrete type so that it is easy to interface them with DSGE models. A `DSGEVAR` type
holds information about the VAR with which we want to combine a given DSGE model and can be easily constructed, given a DSGE object. Once we have constructed a `DSGEVAR` object, then it is straightoforward to estimate the object and compute impulse responses.

```@docs
DSGEVAR
```

### Tips for Using DSGEVARs

* When extensively using DSGE-VARs, we recommend defining your own subspecs in
  `subspecs.jl` because it simplifies the process of saving, reading, and analyzing
  output from estimating and calculating impulse responses for DSGE-VARs.
  See [Advanced Usage](@ref advanced-usage) for a more detailed explanation on changing subspecs.

* The names of the observables must exist as either observables
  or pseudo-observables in the DSGE because for most
  DSGEVAR methods, we need to construct the state space
  representation of the DSGEVAR using information from
  the underlying DSGE.

* It is important to be careful about the order of the
  observables when constructing a DSGEVAR. Whether you define
  the names of the observables by calling `update!` or
  by creating a subspec, we assume that the order of the observables
  corresponds to the observable's index in the data and
  in the state space representation of the DSGEVAR. In the example
  provided above, if we estimate the DSGEVAR on data
  or construct the state space representation of the DSGEVAR,
  we assume that the order of observables in the data array,
  which has dimensions `nobs x nperiods`, is `:obs_gdp`
  in the first row, `:obs_cpi` in the second row, and
  `:obs_nominalrate` in the third row.

* When using functions that use the DSGE as a prior for a VAR (as opposed to
  a VAR approximation of the DSGE), then an intercept term is assumed and cannot
  be turned off. For example, the following two functions computes the VAR coefficients
  and innovations variance-covariance matrix for a `DSGEVAR` object `m`.
  The first one is for a VAR approximation of the DSGE in `m`,
  and it allows the user to specify whether or not they want an intercept term
  using the keyword `use_intercept`. The second function is for using the DSGE
  as a prior for a VAR estimated on `data`. This function does not have the
  `use_intercept` keyword because we require an intercept term when using
  the DSGE as a prior for a VAR.

```
compute_system(m; apply_altpolicy = false,
               regime_switching = false, n_regimes = 2,
               check_system = false, get_system = false,
               get_population_moments = false, use_intercept = false,
               verbose = :high)
compute_system(m, data; apply_altpolicy = false,
               regime_switching = false, n_regimes = 2,
               check_system = false, get_system = false,
               get_population_moments = false,
               verbose = :high)
```
