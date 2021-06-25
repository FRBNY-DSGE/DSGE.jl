# [Forecast Decomposition](@id forecast-decomp)

```@meta
CurrentModule = DSGE
```

Separate from the standard [Forecasting](@ref forecast-step) routines, we have also implemented a function `decompose_forecast` for explaining why a forecast changes as new data becomes available and new estimations are run. Please note that
this function **does not** decompose a forecast into the shocks that produce it. For example, if you want
to understand whether TFP or financial shocks are driving a forecast, then you should
be calculating the shock decomposition output variable (see [Calculating Shock Decompositions](@ref calc-shock-dec)).

```@docs
DSGE.decompose_forecast
```

For an example of how to use this functionality, see [decompose_forecast.jl](https://github.com/FRBNY-DSGE/DSGE.jl/tree/main/examples) on the Github page (or directly inside the directory where DSGE.jl has been downloaded).
