# [Computing Means and Bands](@id means-bands)

```@meta
CurrentModule = DSGE
```

## Procedure

After running a full-distribution forecast, we are often interested in finding
means and density bands of the various forecast outputs. This will allow us to
plot our estimation of the full distribution of the forecast outputs.

**Main Steps:**

- *Load data:* Load and transform data, population data, and population forecast
  (required for forecast output transformations)

- *Read in forecast outputs:* Read in the outputs saved by
  [`forecast_one`](@ref), one variable (e.g. one observable) and one
  `output_type` at a time.

- *Transform forecast outputs:* If necessary, use `reverse_transform` to apply
  transformations specified in the `Observable` or `PseudoObservable` type to
  the given forecast output series.

- *Compute means and bands:* Compute the means and density bands of the forecast
  output.

- *Write to file:* For each `output_var`, Write a `MeansBands` object (see
  [The MeansBands Type](@ref) below) to the file specified by
  `get_meansbands_output_file`.

Computing means and bands is done by calling `compute_meansbands`. If desired, you
can also write your computed means and bands as matrices by calling
`meansbands_matrix_all`.

For example, to compute means and bands for an unconditional, full-distribution
forecast of states and observables:

``` julia
m = AnSchorfheide()
compute_meansbands(m, :mode, :none, [:forecaststates, forecastobs])
```

## Weighted Averages of Full-Distribution Forecasts

Sometimes a forecaster may want to combine several
different "scenarios" to construct a forecast. For example,
there may be uncertainty about the right model to use.
Rather than pick just one model, a forecaster would instead
want to somehow combine the forecasts from all the models.
As another example, a forecaster may not be sure
what policies will be implemented in the future, so
it may be easier to forecast the future by constructing
different plausible scenarios.

We implement forecast combination as a weighted average
of forecast paths drawn from different forecasts.
A full-distribution forecast is just a matrix of different
possible forecast paths, which approximate the true distribution
of forecast paths. To construct a weighted average
of forecasts from different scenarios, we just need to draw randomly from
each scenario's distribution of forecast paths.

For example, suppose `m1`, `m2`, `m3` are three different models
with unconditional forecasts identified by `forecast_string1`, `forecast_string2`,
and `forecast_string3`. Then an equal weighted average of these forecasts
identified by the tag `combo_forecast_string` can be
calculated by running

```
input_types = [:full, :full, :full]
cond_types = [:none, :none, :none]
compute_meansbands([m1, m2, m3], input_types, cond_types, output_vars;
                        weights = [1/3, 1/3, 1/3],
                        forecast_strings = [forecast_string1, forecast_string2, forecast_string3],
                        combo_forecast_string = combo_forecast_string)
```

The results are saved using the file paths implied by the first model in the vector of models
passed as the first input argument (`m1`). To get the resulting `MeansBands`, run

```
mb = read_mb(m1, :full, :none, output_var; forecast_string = combo_forecast_string)
```

## Functions for Calculating Means and Bands

```@docs
compute_meansbands
```

## The `MeansBands` type

``` @docs
MeansBands
```
