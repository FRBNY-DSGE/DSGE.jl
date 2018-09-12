# Macroeconomic Forecasting with DSGEs Using Julia and Parallel Computing

*Marco Del Negro, Abhi Gupta, Pearl Li, Erica Moszkowski*

*April 17, 2017*

In December 2015, we announced [DSGE.jl](https://github.com/FRBNY-DSGE/DSGE.jl),
our open-source, Julia-language package for working with dynamic stochastic
general equilibrium (DSGE) models. At that time, *DSGE.jl* contained only the code
required to specify, solve, and estimate such models using Bayesian
methods. Now, we present the additional code needed to produce economic
forecasts using estimated DSGE models. This new code replicates our MATLAB
codebase while being more efficient, easier to read, and open source.

As we noted in our
[last post](http://libertystreeteconomics.newyorkfed.org/2015/12/the-frbny-dsge-model-meets-julia.html)
and its corresponding
[technical post](https://frbny-dsge.github.io/DSGE.jl/latest/MatlabToJuliaTransition.html),
porting our code to Julia presented us with the opportunity to improve both our
code's performance and our team's workflow. While the estimation step was
largely a direct port, we redesigned the forecast section to obtain code that is
faster and easier to use. In this post, we will discuss the performance
improvements we have achieved in forecasting the DSGE model, as well as the
design principles and Julia tools (particularly related to parallel computing)
that helped us achieve those results.


## Performance Improvements

To motivate our decision to redesign the forecasting code, we first present some
overall performance comparisons between our MATLAB and Julia codebases. Because
the design of the code has changed significantly, these results should *not* be
taken as a horse race between Julia and MATLAB. Rather, they should indicate the
extent to which our design decisions, in conjunction with the power of the Julia
language, have improved the process of running a DSGE forecast.

These tests were conducted on a single core on an Intel® Xeon® E5-2697 v2
2.70GHz CPU running GNU/Linux. The exception is computing all the
full-distribution results, which was done using 50 parallel workers.

**Benchmark Times Relative to MATLAB 2014a (Smaller is Better)**

| Test                                                           | MATLAB (2014a) | Julia (0.4.5) |
| -------------------------------------------------------------- | -------------- | ------------- |
| Simulation smoothing                                           | 1.00           | 0.38          |
| Forecasting                                                    | 1.00           | 0.24          |
| Computing shock decompositions                                 | 1.00           | 0.12          |
| Full set of forecast outputs (modal parameters)                | 1.00           | 0.10          |
| Full set of forecast outputs (full distribution of parameters) | 1.00*          | 0.22          |

*Unlike the other steps being tested, the full-distribution forecast timing was
run in MATLAB 2009a. Our code relies on MATLAB parallelization features that
were deprecated with the introduction of the Parallel Computing Toolbox.

Post estimation, we produce a number of forecast-related outputs, either at the
mode or using the full estimated posterior distribution. The tasks involved
include smoothing, forecasting (both enforcing the zero lower bound and not),
and computing shock decompositions (exercises that allow us to account for the
evolution of observed variables in terms of their driving forces).

With our most recent model, which is available in *DSGE.jl*, we can compute all
the full-distribution forecast outputs in approximately fifteen minutes. In
comparison, the same computations in MATLAB typically take about seventy
minutes. As a result, we can experiment with different options and correct
mistakes much more flexibly than we could previously.

In the next sections, we discuss the design principles that guided our port, as
well as the Julia parallel programming tools that enabled us to write efficient
parallel code.


## Design Principles

Our goal when porting the forecast step to Julia was to write code that could
efficiently produce easy-to-interpret results. Furthermore, because we are often
interested in looking at just one kind of result (for instance, impulse response
functions), we wanted it to be equally simple to produce a single, very specific
output as it was to produce all results. These two goals translated into two
related principles that guided our Julia code development: *type-orientation*
and *modularity*.

**Type-Orientation**

As we discussed in our previous post, Julia's type system allows us to write
very clean, well-structured code, and we use types heavily throughout the
codebase. For example, in the forecast step, we use types heavily to keep track
of the information we need to download and transform our input data. As an
example, let's consider the process of downloading and transforming the GDP
series for use in the DSGE model. First, using the
[FredData.jl](https://github.com/micahjsmith/FredData.jl) package, we pull the
aggregate nominal GDP series (in dollars) from the Federal Reserve Economic
Database ([FRED](https://fred.stlouisfed.org/)) programmatically. Before the
estimation, we transform this series into the appropriate units for the
log-linearized model: quarter-to-quarter log differences of real, per-capita
GDP. After estimating and forecasting the model, we finally transform the
results into the units most frequently discussed by policymakers and
researchers: annualized GDP growth in percentage terms.

We wanted a simple way to keep track of all of the information associated with
the GDP variable in a single place. To do this, we created a new Julia type
called an `Observable`. An instance of the `Observable` type bundles together
the name of the variable, sources used to create the series, and all
transformations associated with that series. An instance of this `Observable`
type has the following fields:

```julia
type Observable
    key::Symbol
    input_series::Vector{Symbol}
    fwd_transform::Function
    rev_transform::Function
    name::UTF8String
    longname::UTF8String
end
```

The `key`, `name`, and `longname` fields serve similar but slightly different
purposes. The `key` is used as the primary way we refer to the GDP variable in
the code: when we construct the entire dataset, we create a `DataFrame`
(2-dimensional table) and label each series with its `key`. By contrast, `name`
is a longer-form name that we intend to use to label plots, while `longname` is
more of a description of the series. This information helps us to label
variables easily and keep the code clear.

The more interesting fields are `input_series`, `fwd_transform`, and
`rev_transform`. The `input_series` field is a vector of `Symbol`s, each of
which must be of the form `:SERIES__SOURCE`. In the case of GDP, this field is
the vector `[:GDP__FRED, :CNP16OV__FRED, :GDPDEF__FRED]`. All of these series
come from FRED, and in particular, we use the nominal GDP, working-age civilian
population, and GDP deflator series to construct the real per-capita GDP growth.

The `fwd_transform` and `rev_transform` fields encode the transformations we
make to the GDP series to go from raw data to model units and from model units
to output units, respectively. These fields are particularly interesting because
they must be populated by objects that are of type `Function`. That's right—a
function is an instance of the `Function` type! Therefore, a given function is
really no different than any other variable in Julia. That means we can define
any function we want (abstract, named, with or without keyword arguments) and
assign the name of that function to the `fwd_transform` and `rev_transform`
fields. In the data step of the code, for instance, we can retrieve the name of
the function by querying the `Observable` object and then apply the function to
an appropriate set of arguments. This is a very direct method of looking up
which transforms to apply, and simultaneously provides the opportunity for us to
abstract common transformations into an appropriately named
function. Abstraction is a technique for encapsulating low-level functionality
or pieces of data into a well-named, reusable function or type. In our case,
abstracting transformations into functions is useful because multiple
observables can make use of the same commonly used functions.

Finally, we can construct the `gdp` observable as follows:

```julia
data_series = [:GDP__FRED, :CNP16OV__FRED, :GDPDEF__FRED]
fwd_transform = function (levels) ... end    # an anonymous function definition
rev_transform = loggrowthtopct_annualized_percapita
obs_gdp = Observable(:obs_gdp, data_series, fwd_transform, rev_transform,
    "Real GDP  Growth", "Real GDP Growth Per Capita")
```

We then store `obs_gdp` in a `Dict{Symbol, Observable}`, a lookup table that
allows us to look up `Observable` objects, which is in turn stored in the
`observables` field of the model object. We can query the model object for the
`rev_transform` of` gdp_obs` by simply calling
`m.observables[:gdp_obs].rev_transform` (where `m` is an instance of a model
type). Since this information is stored inside the model object for every
observable, it is automatically available to every function that accepts a model
object—helping us keep our function calls manageable and our data organized.

We have found Julia's type system to be a helpful way to abstract the details
associated with transforming data to and from various units. `Observable`s are
clearly a *DSGE.jl*-specific example of a user-defined type, but we hope this
discussion illustrates how Julia types and effective abstraction can help
economists structure and clarify their code.

**Modularity**

Most software systems (and economic models, for that matter) are designed to
produce a wide variety of outputs. Macroeconomists often want to produce tables
of parameters, impulse response functions, and time series plots for different
economic variables. Often, users want to choose which of a set of possible
outputs to compute. In a DSGE model, it is common to compute smoothed histories
and forecasts of observables and unobservable states, shock decompositions
(which decompose the path of each economic variable into the shocks responsible
for its fluctuations), and impulse response functions. Additionally, users may
want to change various settings. In our case, we can choose to forecast using
the modal parameters or a selection of draws from the posterior distribution of
the parameters. We can decide whether or not to enforce the zero lower bound on
nominal interest rates. We can use no data from the current quarter, condition
on only financial data from the current quarter, or use both financial data and
GDP data from the current quarter. We can choose from several different
smoothers to compute smoothed histories of states and observables.

Producing and storing all of these results takes both time and disk space. As
users of our own old codebase, we found that these costs were often burdensome
if we only wanted to produce a single result (for instance, an unconditional
shock decomposition for GDP growth). This occurred because the top-level
forecast function always called every subroutine, computed every output, and
returned all outputs. Redesigning the codebase gave us the opportunity to write
code that could produce specific outputs in addition to all outputs.

Fundamentally, in the DSGE forecast, there are three pieces of information we
need to produce the specific outputs desired by the user. First, does the user
want to produce a modal forecast or a full distribution forecast with
uncertainty bands? Second, does she want to condition on any data from the
current quarter? And third, which kinds of outputs does she want to produce
(forecasts, shock decompositions, etc.)?  Once we know the answers to these
questions, we can logically determine which outputs need to be produced and
which can be ignored. Therefore, we present the user with one top-level
function, which takes in these arguments and determines which subroutines need
to be run.

This modular approach to control flow can be taken in any language, but it is an
important component of developing a large software system or economic model and
thus we decided it was important to mention. Writing modular, type-oriented, and
well-abstracted code improves the robustness of our workflow by making our code
and results easier to interpret and less prone to error. In the next section,
we'll discuss the main reason our Julia codebase is so fast: we are able to
exploit Julia's parallel programming tools.


## Parallel Computing

The types of forecast-related computations we do are naturally suited to
parallelization. While our MATLAB code was parallelized to an extent (and was
written before the advent of the MATLAB Parallel Computing Toolbox!), we decided
to reassess our design when we ported the forecast step to Julia. We considered
two approaches: "parallel maps" and distributed storage. The first is largely
similar to our MATLAB parallelization implementation, while the latter takes
advantage of the
[DistributedArrays.jl](https://github.com/JuliaParallel/DistributedArrays.jl)
package and represents a substantial design shift. Over the course of
development, we learned a great deal about writing effective parallelized Julia
code and about parallel computing in general. Though the distributed storage
approach did not end up improving on the parallel mapping approach, our final
Julia code is faster and better designed than the original MATLAB
implementation.

Like many academic institutions, the New York Fed's Research Group maintains a
Linux-based cluster for use by the economists and RAs. This setup allows us to
distribute computing jobs across multiple processes on multiple compute nodes,
so that non-serially dependent jobs can be executed at the same time. However,
our jobs must also coexist with those of other researchers, which limits both
the amount of CPU time and memory we can use before disrupting other work. Our
code is designed to take advantage of the features of and respect the
constraints of this environment.

During the estimation step, we simulate drawing a large number of parameters
(typically 100,000) from their posterior distribution. In the forecast step,
these draws are read in and used to compute the desired outputs for our observed
variables and the latent states. As discussed before, these outputs can include
smoothed shock times series, forecasts, shock decompositions, and impulse
response functions. Since these computations are independent for each parameter
draw, forecasting using the full distribution lends itself well to
parallelization.

To reduce our impact on other users of the cluster, we make use of a "blocking"
scheme in our Julia code. The parameter draws are read in blocks of typically
5,000 draws. These draws are then immediately distributed using Julia's `pmap`
function ("parallel map") to worker processes, each of which carries out the
entire forecast step for just one draw. When all of the draws from that block
have completed, the originator process re-collects the forecast outputs and
saves them to disk. This repeats until all blocks have completed. Through this
blocking, we can avoid keeping too much data in memory, since we only operate on
a fraction of the parameter draws at any given time. However, we can still write
structured output files using the HDF5 file format, which allows us to write to
specific subsets of pre-allocated arrays, so that the end result is as if we had
computed all the draws at once without blocking.

Before settling on this version, we also tried using Julia's
[DistributedArrays.jl](https://github.com/JuliaParallel/DistributedArrays.jl)
package, which distributes large arrays in memory over many processes. This
allowed us to hold all of our parameter draws and their corresponding forecast
outputs in memory at the same time, and it allowed each process to operate on
the parameter draws it held locally without needing to copy data back and forth
between processes. However, using distributed arrays also forced us to
explicitly handle lower-level tasks like assigning parameter draws to
processes. Since each process handled a predetermined set of draws, it was not
easy to reallocate draws if some of the compute nodes on which the processes
lived happened to be busier than others on a particular day. Switching to `pmap`
allowed us to abstract away from many of these concerns, as it has already been
optimized to take advantage of the aforementioned independence of parameter
draws.


## StateSpaceRoutines.jl

A big benefit of using Julia is the large and growing package ecosystem, which
allows all users to access high-quality open-source code. Thanks to this system,
Julia developers can focus their development time on the issues and projects
they really care about, without having to repeatedly reinvent the wheel. For
example, *DSGE.jl* depends on the
[DataFrames.jl](https://github.com/JuliaStats/DataFrames.jl) package to help us
manage data and dates. Similarly, *DSGE.jl* is available for members of the
community to modify, extend, and make use of as they see fit. In this spirit, we
have decided to break out some DSGE-independent components of *DSGE.jl* into
their own package.

DSGE models define a linear system that links observed variables to unobserved
states. In order to actually perform inference on these latent states, we apply
the Kalman filter and smoothing algorithms. State space models are commonly used
across many disciplines, and indeed the routines we use in *DSGE.jl* can be
applied to any sort of linear state space model. As such, we have decided to
move the filtering and smoothing routines that we have historically used with
the DSGE model into
[StateSpaceRoutines.jl](https://github.com/FRBNY-DSGE/StateSpaceRoutines.jl), a
new package that will provide *DSGE.jl*-independent filtering and smoothing
routines.

*StateSpaceRoutines.jl* currently features the one filter and four smoothers we
most commonly use in *DSGE.jl*. On the filtering front, we implement the
standard Kalman filter found in James Hamilton's *Time Series Analysis*
*StateSpaceRoutines.jl* also contains two Kalman smoothers and two simulation
smoothers. In addition to the Kalman smoother presented in *Time Series
Analysis*, we also have Jan Koopman's disturbance smoother from his paper
*Disturbance Smoother for State Space Models*. The two simulation smoothers are
based on Carter and Kohn's *On Gibbs Sampling for State Space Models* and Durbin
and Koopman's *A Simple and Efficient Simulation Smoother for State Space Time
Series Analysis*. In our experience, the Koopman smoother is faster than the
standard Kalman smoother, as it does not require us to calculate the pseudo-
inverses of the predicted variance matrices. For the same reason, we have also
found that the Durbin and Koopman simulation smoother is faster than the Carter
and Kohn one. All of these methods support time-varying matrices and
variances. We use this feature to model pre–zero-lower-bound and
zero-lower-bound regimes in our DSGE models, but the functionality is general
enough to be applied to a wider range of models with regime switching or time
varying matrices and variances. We hope that the broader Julia community finds
these functions as useful as we have!


## Disclaimer

This post reflects the experience of the authors with Julia and MATLAB and does
not represent an endorsement by the Federal Reserve Bank of New York or the
Federal Reserve System of any particular product or service. The views
expressed in this post are those of the authors and do not necessarily reflect
the position of the Federal Reserve Bank of New York or the Federal Reserve
System. Any errors or omissions are the responsibility of the authors.


## References

Carter, C. and Cohn, R. (1994). On Gibbs Sampling for State Space models. *Biometrika*.

Durbin, K. and Koopman, S. (2002). A Simple and Efficient Smoother for State Space Time Series Analysis. *Biometrika*.

Hamilton, J. (1994). *Time Series Analysis*. Princeton: Princeton University Press.

Koopman, S. (1993). Disturbance Smoother for State Space Models. *Biometrika*.
