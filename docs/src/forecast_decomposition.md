# [Forecast Decomposition](@id forecast-decomp)

```@meta
CurrentModule = DSGE
```

Separate from the standard [Forecasting](@ref forecast-step) routines, we have also implemented a function `decompose_forecast` for explaining why a forecast changes as new data becomes available and new estimations are run.

```@docs
DSGE.decompose_forecast
```

For an example of how to use this functionality, see [decompose_forecast.jl](https://github.com/FRBNY-DSGE/DSGE.jl/tree/main/examples) on the Github page (or directly inside the directory where DSGE.jl has been downloaded).
