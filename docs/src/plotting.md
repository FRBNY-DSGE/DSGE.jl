# [Plotting](@id plotting-doc)

``` @meta
CurrentModule = DSGE
```
The DSGE plotting code uses [Plots](https://juliaplots.github.io/). We
typically use GR and Plotly as
[backends](https://juliaplots.github.io/backends/); however, the goal of Plots
is that all backends should be supported interchangeably. In each of the functions
listed below, there are methods which take in an `AbstractModel` and those which
take in some lower-level input arguments (typically one or more
`MeansBands`). See the individual function docstring for details.

## Plotting Estimation Results

- `plot_prior_posterior`: plot the prior distribution overlaid on a histogram of
  posterior draws

## Plotting Forecasts

- `plot_history_and_forecast`: plot a historical and forecasted series, possibly
  with uncertainty bands (if for a full-distribution forecast)
- `plot_forecast_comparison`: plot two sets of histories and forecasts in one plot
- `hair_plot`: plot many forecasts as "hairs" coming out of some realized data
  series
- `plot_shock_decomposition`: plot the contributions of individual shocks as
  a bar plot, with a line for the detrended mean forecast
- `plot_impulse_response`: plot impulse response functions

## Other Plots

- `plot_altpolicies`: plot forecasts under several alternative policies in one
  plot
- `plot_scenario`: plot a forecast conditional on some alternative scenario, in
  deviations from some baseline
