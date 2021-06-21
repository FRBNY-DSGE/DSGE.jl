# Learning How to Use DSGE.jl

DSGE.jl is designed to facilitate the fast and efficient estimation and forecasting of DSGE models.
For this reason, DSGE.jl does not have as friendly an API as other packages like Dynare.
For example, building new models in DSGE.jl will be more involved than the model script approach
taken by Dynare, and you will need to know basic Julia first.
While the FRBNY DSGE team tries it best to write good documentation,
users may still need to look at source code to learn how some functions work
and to diagnose errors. Please file an [issue](https://github.com/FRBNY-DSGE/DSGE.jl/issues)
or contact one of the developers (see [Project.toml](https://github.com/FRBNY-DSGE/DSGE.jl/tree/main/Project.toml))
if you ever have questions about the package or encounter errors that you are not sure how to fix.
However, try to thoroughly examine the documentation first because, in our experience,
many problems users have had were easily answered by referring them to the correct section of the documentation.

To get started, we recommend taking a look at our
[example scripts](https://github.com/FRBNY-DSGE/DSGE.jl/tree/main/examples)
as you read the documentation because the examples will make it easier to comprehend the documentation.
The two scripts that will be most useful for learning the package are
- [`run_default.jl`](https://github.com/FRBNY-DSGE/DSGE.jl/tree/main/examples/run_default.jl): a simple example of the standard workflow with DSGE.jl
- [`make_packet.jl`](https://github.com/FRBNY-DSGE/DSGE.jl/tree/main/examples/make_packet.jl): auto-generate a packet of plots and figures which
  help the user analyze estimation, forecast, and impulse response results. This script
  also provides an example of how we recommend structuring "main" files that launch
  a forecast and generate results with one "click."
For details on all the example scripts,
see [Running Existing Models](https://frbny-dsge.github.io/DSGE.jl/stable/running_existing_model/).

Below, we give some additional pointers depending on what the user wants from the package.

## Building New DSGE Models

For beginners, we recommend looking at the
[source code](https://github.com/FRBNY-DSGE/DSGE.jl/tree/main/src/models/representative/an_schorfheide/)
for `AnSchorfheide` to see the most basic example of a DSGE model in DSGE.jl.
Unless you are experienced with object-oriented programming and implementing new types in Julia,
we *highly recommend* you start with `AnSchorfheide`. In our experience, the typical economist
who wants to use our code encounters a substantial learning curve, especially if they
try to build their own DSGE models by modifying `Model1002`, the current New York Fed DSGE model.
Because `Model1002` continues to undergo development by the FRBNY DSGE team, it includes
numerous features that require deep knowledge of the package and is generally *not* suitable for
learning to build a new model.

To help users, the source code for `AnSchorfheide` should have additional comments indicating which
functions and lines should be edited when a user wants to build a new model.
For this reason, we recommend that the user copies the source code for the model into a new folder
and edit that source code to develop their new model. Substantial portions of the code
for the model objects in DSGE.jl can be recycled for most DSGE models. Copy and paste
will minimize errors and also reduce the time it takes to get a new model implemented.

For more advanced programmers and for users familiar with the implementation of `AnSchorfheide`,
`SmetsWouters` or `SmetsWoutersOrig` are suitable examples for building medium-scale DSGEs.
The former uses our preferred parameterization of the Smets and Wouters model, and the latter
uses the parameterization in the original paper. A more advanced model that is
closer to the New York Fed DSGE model is `Model904`.
To see the source code of these models, look in the
[representative agent DSGE folder](https://github.com/FRBNY-DSGE/DSGE.jl/tree/main/src/models/representative/).
We recommend examining these models first
before looking at `Model1002` because there are notable differences in the size of these models
compared to `AnSchorfheide`.

Finally, after learning our implementation of medium-scale DSGEs, the user may start to look at
`Model1002` if they want to base their new model off of the New York Fed DSGE model.
In addition, the most advanced features of DSGE.jl are typically implemented only for `Model1002`.
To learn these features, you will likely need to examine the
[source code for `Model1002`](https://github.com/FRBNY-DSGE/DSGE.jl/tree/main/src/models/representative/m1002/).
For example, we have extended `Model1002` to model the economic impact of COVID-19
and the Federal Reserve's policy switch to flexible average inflation targeting. None of the other
DSGE models in DSGE.jl has been extended to include these features, although it is certainly possible to do so.

The user will also want to look at [Solving the Model](@ref solving-dsge-doc) because it will explain
the format in which equilibrium conditions must be provided in order for a DSGE model to be solvable by DSGE.jl.
In particular, this page will guide the user when writing their own model's `eqcond.jl` script.

## Estimation

The user should make sure they have a FRED API key (see these [instructions](https://github.com/micahjsmith/FredData.jl))
if they want to estimate a DSGE model because we load most US data automatically by querying FRED.
Our package is designed to make the data retrieval process as automated as possible, but
this approach means that more setup is required to load data into Julia than simply creating your own
csv file and reading it into memory. Please see [Input Data](@ref input-data-step) for guidance.

A relatively high-level overview of the estimation process is given by
[Estimation](@ref estimation-step). The [`run_default.jl`](https://github.com/FRBNY-DSGE/DSGE.jl/tree/main/examples/run_default.jl)
script will be a particularly useful
example, followed by `test_smc.jl` if the user wants to also learn how to call the Sequential Monte Carlo
sampling algorithm instead of Random-Walk Metropolis-Hastings. The default settings
for estimation are designed for the FRBNY DSGE team's typical work, but most users will
need to adjust the defaults. For guidance on what settings can be changed and what they do,
see [Working with Settings](@ref working-with-settings) as well as
[defaults.jl](https://github.com/FRBNY-DSGE/DSGE.jl/blob/main/src/defaults.jl), which
holds the complete description of default settings for every DSGE model in DSGE.jl. Note, however,
that some of the values of these default settings will be over-written during instantiation of a DSGE object
(see the `model_settings!` function for `AnSchorfheide` in
[`an_schorfheide.jl`](https://github.com/FRBNY-DSGE/DSGE.jl/blob/main/src/models/representative/an_schorfheide.jl)).

## Forecasting

The user should make sure they have a FRED API key (see these [instructions](https://github.com/micahjsmith/FredData.jl))
if they want to forecast a DSGE model because we load most US data automatically by querying FRED.
Our package is designed to make the data retrieval process as automated as possible, but
this approach means that more setup is required to load data into Julia than simply creating your own
csv file and reading it into memory. Please see [Input Data](@ref input-data-step) for guidance.

A relatively high-level overview of forecasting process is given by
[Forecasting](@ref forecast-step) and [Computing Means and Bands](@ref means-bands).
The first section describes the actual "running" of forecasts while the second section
describes the `MeansBands` object we typically utilize to examine forecast output. In particular,
most use cases will not require access to the matrix describing every single sample from
the posterior forecast distribution. Instead, a user only wants to know what the mean forecast
and the uncertainty bands around that mean are. For this reason, the `MeansBands` object
is typically the object users will want in order to analyze the output of their forecast.
A good example script of forecasting is
[`make_packet.jl`](https://github.com/FRBNY-DSGE/DSGE.jl/tree/main/examples/make_packet.jl),
which provides an instance of a typical workflow when forecasting a DSGE model.


For more details on specific forecast outputs that can be computed, the user will find
[Impulse Response Functions](@ref irf-doc), [Alternative Policies](@ref altpol-doc),
[Forecasting Decomposition](@ref forecast-decomp), and [Plotting](@ref plotting-doc) to be informative sections.
Once the user is familiar with the basic applications of the forecasting tools in DSGE.jl,
they may want to look at [Advanced Usage](@ref advanced-usage) for the many special features we have
implemented for forecasting DSGEs.