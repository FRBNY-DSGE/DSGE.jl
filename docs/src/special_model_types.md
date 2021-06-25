# Special Model Types

```@meta
CurrentModule = DSGE
```

In addition to the `AbstractDSGEModel` type, DSGE.jl has several special types that implement
models which are designed to interface with DSGEs. New users should skip this section and return later.

## The `PoolModel` Type
Unlike the other models contained in DSGE, the `PoolModel` type is not a proper DSGE model.
It is a wrapper object for different methods to perform model averaging for two different models,
which do not have to be DSGE models.
For example, a user could average two different vector auto-regressions.
Generally, a user only needs to provide the predictive density scores
of the two models that the user wants to average. The reason is that we
treat the predictive density scores as non-FRED observables. This approach
makes interfacing with the rest of the machinery provided by DSGE.jl very simple.

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

For more details on the theory and performance, see
[Del Negro and Schorfheide (2004)](https://onlinelibrary.wiley.com/doi/full/10.1111/j.1468-2354.2004.00139.x),
[Del Negro and Schorfheide (2006)](https://www.newyorkfed.org/medialibrary/media/research/economists/delnegro/erq206_delnegro.pdf),
[Del Negro, Schorfheide, Smets, and Wouters (2007)][https://www.jstor.org/stable/27638915], and
[Del Negro and Schorfheide (2009)](https://www.aeaweb.org/articles?id=10.1257/aer.99.4.1415).

We implement DSGE-VARs with the `DSGEVAR` concrete type so that it is easy to interface them with DSGE models. A `DSGEVAR` type
holds information about the VAR with which we want to combine a given DSGE model and can be easily constructed, given a DSGE object. Once we have constructed a `DSGEVAR` object, then it is straightforward to estimate the object and compute impulse responses.

```@docs
DSGEVAR
```

### [Tips for Using `DSGEVAR`](@id tips-dsgevar)

* When extensively using DSGE-VARs, we recommend defining your own subspecs in
  `subspecs.jl` because it simplifies the process of saving, reading, and analyzing
  output from estimating and calculating impulse responses for DSGE-VARs.
  See [Advanced Usage](@ref advanced-usage) for a more detailed explanation on changing subspecs.

* The names of the observables must exist as either observables
  or pseudo-observables in the DSGE because for most
  `DSGEVAR` methods, we need to construct the state space
  representation of the `DSGEVAR` using information from
  the underlying DSGE.

* It is important to be careful about the order of the
  observables when constructing a `DSGEVAR`. Whether you define
  the names of the observables by calling `update!` or
  by creating a subspec, we assume that the order of the observables
  corresponds to the observable's index in the data and
  in the state space representation of the `DSGEVAR`. In the example
  provided above, if we estimate the `DSGEVAR` on data
  or construct the state space representation of the `DSGEVAR`,
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

## DSGE-VECMs and the `DSGEVECM` Type

We can extend DSGE-VARs to permit cointegrating relationships between observables
using DSGE-VECMs. A VECM is a [vector error-correction model](https://en.wikipedia.org/wiki/Error_correction_model),
which extend VARs to account for long-run stochastic trends, i.e. cointegration..

For more details on the theory and performance, see
[Del Negro and Schorfheide (2004)](https://onlinelibrary.wiley.com/doi/full/10.1111/j.1468-2354.2004.00139.x),
[Del Negro and Schorfheide (2006)](https://www.newyorkfed.org/medialibrary/media/research/economists/delnegro/erq206_delnegro.pdf), and
[Del Negro, Schorfheide, Smets, and Wouters (2007)][https://www.jstor.org/stable/27638915].

We implement DSGE-VECMs with the `DSGEVECM` concrete type so that it is easy to interface them with DSGE models. The `DSGEVECM`
has very similar behavior to the `DSGEVAR` type, with extensions as needed. For example, the `DSGEVECM` type includes additional fields
to hold information about cointegrating relationships.

```@docs
DSGEVECM
```

### Tips for Using `DSGEVECM`

* The same [tips for `DSGEVAR` models](@ref tips-dsgevar) generally apply for `DSGEVECM` models.

* The names of cointegrating relationships in the field `cointegrating`
  must exist as either observables or pseudo-observables in the DSGE. The reason is
  the same as the reason for why observables must be held.

* In the state space representation of the underlying DSGE corresponding to
  a `DSGE-VECM`, cointegrating relationships are ordered after observables. For example,
  consider the measurement matrix `ZZ`. The first `n_observables` rows correspond to
  the `observables` in `DSGE-VECM`, and the next
  `n_observables + 1:n_cointegrating + n_observables` rows correspond to
  `cointegrating` in `DSGE-VECM`.

* When calculating the `VECM` coefficients of a DSGE-VECM,
  the coefficients are ordered with cointegrating relationships first, followed by
  the intercept term, and concluding with lags of past differences. See the
  docstring of `vecm_approx_state_space`.

* Some cointegrating relationships do not need to be added to the measurement matrix
  in the state space representation of a DSGE model. These relationships are considered
  "additional" ones and are added to the field `cointegrating_add`. To compute the constant
  vector which specify these additional relationships, we use `compute_DD_coint_add`.
  See its docstring for notes on usage.

## Auxiliary Methods for DSGE-VARs and DSGE-VECMs
Listed below are some methods used internally in but not exported by DSGE.jl
that users may also find useful. We also recommend looking at the various utility
functions in `abstractvarmodel.jl`. Many of these functions are wrappers for similarly
named functions defined on `AbstractDSGEModel` objects.

```@autodocs
Modules = [DSGE]
Pages = ["models/var/util.jl"]
Order = [:function]
```
```@docs
var_approx_state_space
vecm_approx_state_space
compute_DD_coint_add
measurement_error
dsgevar_likelihood
dsgevecm_likelihood
```
