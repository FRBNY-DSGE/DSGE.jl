# [Impulse Responses](@id irf-doc)

```@meta
CurrentModule = DSGE
```

We provide many different types of impulse responses for DSGEs,
VARs, VECMs, DSGE-VARs, and DSGE-VECMs. The forecast step allows the user to automatically compute
"structural" impulse responses specifically for DSGEs,
but for some purposes, a user may just want impulse responses
without having to compute any other quantities. We provide this functionality
with the `impulse_responses` function. See the end of this page for the docstrings
of all available impulse response functions.

We overload impulse_responses to cover specific use cases. For any
`AbstractDSGEModel`, we can compute the impulse responses over
a specified horizon
for all endogenous state variables, observables, and pseudo-observables
by running

```julia
m = AnSchorfheide()
system = compute_system(m)
horizon = 40
states_irf, obs_irf, pseudo_irf = impulse_response(system, horizon)
```

For an `AbstractRepModel` (a subtype of `AbstractDSGEModel` for representative
agent models), we can also grab the impulse responses by running


```julia
states_irf, obs_irf, pseudo_irf = impulse_response(m, system)
```

This use case requires the user to add a setting
under the key `:impulse_response_horizons`, which is
set by default to 40.


If a user wants to specify a subset of the exogenous shocks
and the size of those shocks, the user can run

```julia
shock_names = [:g_sh, :b_sh]
shock_values = [1.0, 1.0]
impulse_responses(m, system, horizon, shock_names, shock_values)
```

For the response of an endogenous state or observable to a specific shock,

```julia
shock_name  =  :g_sh
var_name = :obs_gdp
var_value = 0.
impulse_responses(m, system, horizon, shock_name , var_name, var_value)
```

## DSGE Impulse Responses
There are two categories of impulse responses for DSGEs provided by DSGE.jl.
It is easy to distinguish them by examining the state space form of a DSGE model (see [Solving](@ref solving-dsge-doc)):
```math
\begin{aligned}
s_t &= T s_{t-1} + R \epsilon_t + C & \epsilon_t &\sim N(0, Q) & \mathrm{(transition)} \\
y_t &= Z s_t + D + u_t & u_t &\sim N(0, E) & \mathrm{(measurement)}
\end{aligned}
```
Impulse responses in the first category are "structural" impulse responses, which are
the response of states and observables to the exogenous structural shocks ``\epsilon_t``.


Impulse responses in the second category are "observables-identified" impulse responses.
First, we may suppose that the measurement equation generically follows

```math
y_t = F(s_t) + \eta_t,
```

where ``F(\cdot)`` is some function of the unknown states ``s_t``,
and ``\eta_t`` are random innovation to the observables ``y_t``.
By innovations, we mean that these random variables are potentially
endogenous shocks, i.e. shocks which do not have a causal interpretation.
An "observables-identified" impulse response specifies
a certain response of ``y_t`` to the innovations ``\eta_t``,
and uses this response to identify the underlying structural shocks which are
consistent with these innovations. A DSGE identifies these innovations
using the state space form of a DSGE.

We provide three types of "observables-identified" impulse responses for DSGEs.

- Short-Run Cholesky
- Long-Run Cholesky
- Maximizing Explained Cyclical Variance

We document the details of the identification in the docstrings of these impulse
response functions. For the first two types of impulse responses,
search the docstrings at the end of the page for

```
function impulse_responses(system::System{S}, horizon::Int, permute_mat::AbstractMatrix{T},
                           shocks::AbstractVector{S} = Vector{S}(undef, 0);
                           restriction::Symbol = :short_run, flip_shocks::Bool = false,
                           get_shocks::Bool = false) where {S <: Real, T <: Number}
```

For the third type of impulse response, search for

```
function impulse_responses(system::System{S}, horizon::Int, frequency_band::Tuple{S,S},
                           n_obs_shock::Int; flip_shocks::Bool = false,
                           get_shocks::Bool = false) where {S <: Real}
```


## VAR Impulse Responses
While we have not yet implemented a VAR model, we do have impulse
responses often used on VARs because of DSGE-VARs. Consider the VAR

```math
y_t = X_t \beta + \epsilon_t,
```

where ``X_t`` is a matrix of the lags of ``y_t``, ``\beta`` are the
VAR coefficients, and ``\epsilon_t \sim N(0, \Omega)`` are the innovations
to observables.

We provide three types of impulse responses, each of which
provide a different way of identifying orthogonalized shocks
from the innovations.

- Short-Run Cholesky
- Long-Run Cholesky
- Maximizing Explained Cyclical Variance

These impulse responses are named the same as the observables-identified
impulse responses for DSGEs because they are considering the same
response of observables to the innovations.
However, the treatment of identification is different when using a VAR
because the mathematical structure of a VAR is not the same as a DSGE's.
As a result, there are slight differences between these impulse responses
and the observables-identified DSGE impulse responses.


```
impulse_responses(β, Σ, n_obs_shock, horizon, shock_size = 1;
    method = :cholesky, flip_shocks = false, use_intercept = true,
    frequency_band = (2π/32, 2π/6)) where {S<:Real}
```



## DSGE-VAR Impulse Responses

There are two types of impulse responses we can compute for a DSGE-VAR.
For both types, we first draw from the posterior distributions of the ``VAR(p)`` coefficients
and innovations variance-covariance matrix, where ``p`` is
the number of lags. With these draws, we can do one of two things:

1. Compute the VAR impulse response implied by the draws.
2. Use the DSGE's structural impact response (i.e. the first period of an impulse response)
   to identify a mapping from the (endogenous) innovations in the VAR
   to the structural shocks of the DSGE.

The first type of impulse response uses the same code as the VAR impulse responses
once we compute the coefficients and innovations variance-covariance matrix.
We call the second type of impulse responses
"DSGE-VAR rotation impulse responses" because we effectively use the DSGE
to identify a rotation matrix.

For the first type of impulse response, see

```
function impulse_responses(m::AbstractDSGEVARModel{S}, data::AbstractArray{S}, method::Symbol,
                           n_obs_shock::Int; horizon::Int = 0 ,use_intercept::Bool = false,
                           flip_shocks::Bool = false, verbose::Symbol = :none) where {S <: Real}

function impulse_responses(m::AbstractDSGEVARModel{S}, method::Symbol,
                           n_obs_shock::Int; horizon::Int = 0 ,use_intercept::Bool = false,
                           flip_shocks::Bool = false, verbose::Symbol = :none) where {S <: Real}
```

The second function is for the specific case when ``\lambda = \infty``, where the data does not matter for the impulse response.

For the second type of impulse responses, see

```
function impulse_responses(m::AbstractDSGEVARModel{S}, data::AbstractArray{S},
    X̂::Matrix{S} = Matrix{S}(undef, 0, 0);
    horizon::Int = 0, MM::Matrix{S} = Matrix{S}(undef, 0, 0),
    flip_shocks::Bool = false, draw_shocks::Bool = false,
    deviations::Bool = false,
    verbose::Symbol = :none) where {S <: Real}
```


## VECM Impulse Responses
While we have not yet implemented a VECM model, we do have impulse
responses often used on VECMs because of DSGE-VECMs. Consider the VECM

```math
\Delta y_t = e_{t - 1} \beta_e +  X_t \beta_v + \epsilon_t,
```

where ``e_{t - 1}`` are cointegrating relationships (the error correction terms);
``X_t`` is a matrix of the lags of ``\Delta y_t``; ``\beta_e`` and ``\beta_v`` are the
VECM coefficients; and ``\epsilon_t \sim N(0, \Omega)`` are the innovations
to observables. We identify orthogonalized shocks for VECMs
from the innovations using the short-run Cholesky method. Other methods
have yet to be implemented, hence passing keywords specific to these methods
will not do anything (namely `frequency_band`).
Find the docstring of the following function for details.

```
impulse_responses(β, Σ, coint_mat, n_obs_shock, horizon, shock_size = 1;
    method = :cholesky, flip_shocks = false, use_intercept = true,
    frequency_band = (2π/32, 2π/6)) where {S<:Real}
```


## DSGE-VECM Impulse Responses

Most of the impulse responses for DSGE-VARs have been implemented for DSGE-VECMs.
The two impulse responses that have not been implemented, due to the differences in
VECMs and VARs, are the `maxBC` and `cholesky_long_run` impulse responses.
For the first type of impulse responses, which use the VECM impulse response code,
see the functions

```
function impulse_responses(m::AbstractDSGEVECMModel{S}, data::AbstractArray{S},
                           coint_mat::AbstractMatrix{S}, method::Symbol,
                           n_obs_shock::Int; horizon::Int = 0 ,use_intercept::Bool = false,
                           flip_shocks::Bool = false, verbose::Symbol = :none) where {S <: Real}

function impulse_responses(m::AbstractDSGEVECMModel{S}, coint_mat::AbstractMatrix{S}, method::Symbol,
                           n_obs_shock::Int; horizon::Int = 0 ,use_intercept::Bool = false,
                           flip_shocks::Bool = false, verbose::Symbol = :none) where {S <: Real}
```

The second function is for the specific case when ``\lambda = \infty``, where the data does not matter for the impulse response.

For the second type of impulse responses, which we call rotation impulse responses, see

```
function impulse_responses(m::AbstractDSGEVECMModel{S}, data::AbstractArray{S},
                           coint_mat::AbstractMatrix{S},
                           X̂::Matrix{S} = Matrix{S}(undef, 0, 0);
                           horizon::Int = 0, MM::Matrix{S} = Matrix{S}(undef, 0, 0),
                           flip_shocks::Bool = false, draw_shocks::Bool = false,
                           deviations::Bool = false,
                           verbose::Symbol = :none) where {S <: Real}
```

## Wrappers for Impulse Response Functions

The `forecast_one` function provides a wrapper for computing
structural DSGE impulse responses when drawing from a distribution of parameters
and for saving these impulse responses as `MeansBands` objects
(see [Computing Means and Bands](@ref means-bands)).

However, `forecast_one` currently cannot
compute observables-identified DSGE impulse responses, VAR impulse responses,
or DSGE-VAR impulse responses, and we do not plan on modifying `forecast_one`
to make it possible to do so. Instead, we provide three wrapper functions
specifically for computing means and bands for these impulse responses.
Please see their docstrings for details. The first wrapper is for
observables-identified DSGE impulse responses, and the second two
are for DSGE-VAR impulse responses. The first one applies generically to
a DSGE-VAR with any prior weight ``\lambda``, but the second one
is a convenience wrapper for the case of ``\lambda = \infty``,
which is equivalent to computing the impulse responses of the
VAR approximation to a DSGE.

No wrappers for DSGE-VECM impulse responses have been implemented
because we have not constructed a DSGE model that can be
interfaced with the `DSGEVECM` type yet. As a result, there are
no explicit test cases for these wrappers. We have decided
against implementing wrappers for DSGE-VECM impulse responses until
we have test cases to guarantee the wrappers do not have bugs.

For observables-identified DSGE impulse responses, find

```
function impulse_responses(m::AbstractDSGEModel, paras::Matrix{S},
                           input_type::Symbol, method::Symbol, n_obs_shock::Int,
                           output_vars::Vector{Symbol} =
                           [:irfstates, :irfobs, :irfpseudo]; parallel::Bool = false,
                           permute_mat::Matrix{S} = Matrix{Float64}(undef,0,0),
                           frequency_band::Tuple{S,S} = (2*π/32, 2*π/6),
                           flip_shocks::Bool = false,
                           density_bands::Vector{Float64} = [.5, .6, .7, .8, .9],
                           create_meansbands::Bool = false, test_meansbands::Bool = false,
                           minimize::Bool = true,
                           forecast_string::String = "",
                           do_rev_transform::Bool = false,
                           verbose::Symbol = :high) where {S<:Real}
```

For DSGE-VAR impulse responses, find

```
function impulse_responses(m::AbstractDSGEVARModel{S}, paras::Matrix{S},
                           data::Matrix{S}, input_type::Symbol, method::Symbol;
                           parallel::Bool = false,
                           frequency_band::Tuple{S,S} = (2*π/32, 2*π/6),
                           n_obs_shock::Int = 1, draw_shocks::Bool = false,
                           flip_shocks::Bool = false,
                           X̂::AbstractMatrix{S} = Matrix{S}(undef, 0, 0),
                           deviations::Bool = false
                           density_bands::Vector{Float64} = [.5, .6, .7, .8, .9],
                           create_meansbands::Bool = false, test_meansbands::Bool = false,
                           minimize::Bool = true,
                           forecast_string::String = "",
                           verbose::Symbol = :high) where {S<:Real}
```


For the impulse response of a VAR approximation to a DSGE, find

```
function impulse_responses(m::AbstractDSGEModel, paras::Union{Vector{S}, Matrix{S}},
                           input_type::Symbol, method::Symbol,
                           lags::Int, observables::Vector{Symbol},
                           shocks::Vector{Symbol},
                           n_obs_shock::Int; parallel::Bool = false,
                           frequency_band::Tuple{S,S} = (2*π/32, 2*π/6),
                           flip_shocks::Bool = false,
                           use_intercept::Bool = false,
                           density_bands::Vector{Float64} = [.5, .6, .7, .8, .9],
                           create_meansbands::Bool = false,
                           minimize::Bool = true,
                           forecast_string::String = "",
                           verbose::Symbol = :high) where {S<:Real}
```

## Docstrings

```@docs
DSGE.impulse_responses
```
