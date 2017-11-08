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

```@docs
compute_meansbands
```


## The `MeansBands` Type

``` @docs
MeansBands
```
