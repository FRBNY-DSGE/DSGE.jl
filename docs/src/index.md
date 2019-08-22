# DSGE.jl

```@meta
CurrentModule = DSGE
```

The *DSGE.jl* package implements the New York Fed DSGE model and provides
general code to estimate many user-specified DSGE models. The package is
introduced in the Liberty Street Economics blog post
[The FRBNY DSGE Model Meets Julia](http://libertystreeteconomics.newyorkfed.org/2015/12/the-frbny-dsge-model-meets-julia.html).

This Julia-language implementation mirrors the MATLAB code included in the
Liberty Street Economics blog post
[The FRBNY DSGE Model Forecast](http://libertystreeteconomics.newyorkfed.org/2015/05/the-frbny-dsge-model-forecast-april-2015.html).

The New York Fed DSGE team is currently working on adding methods to solve nonlinear and
heterogeneous agent DSGE models. Extensions of the DSGE model code may be released
in the future at the discretion of the New York Fed.

## Table of Contents

```@contents
Pages = [
  "model_design.md",
  "running_existing_model.md",
  "advanced_usage.md",
  "input_data.md",
  "frbny_data.md",
  "implementation_details.md",
  "solving.md",
  "estimation.md",
  "forecast.md",
  "means_bands.md",
  "altpolicy.md",
  "scenarios.md",
  "plotting.md",
  "algorithms.md",
  "contributing.md",
  "MatlabToJuliaTransition.md",
  "julia_forecasting.md",
  "license.md"
]
```

## Acknowledgments

Developers of this package at the
[New York Fed](https://www.newyorkfed.org/research) include

* [Michael Cai](https://github.com/caimichael)
* [William Chen](https://github.com/chenwilliam77)
* [Abhi Gupta](https://github.com/abhig94)
* [Pearl Li](https://github.com/pearlzli)
* [Ethan Matlin](https://github.com/ethanmatlin)
* [Erica Moszkowski](https://github.com/emoszkowski)
* [Reca Sarfati](https://github.com/rsarfati)
* [Micah Smith](https://github.com/micahjsmith)

Contributors to this package at [QuantEcon](http://quantecon.org) include

* [Zac Cranko](https://github.com/ZacCranko)
* [Spencer Lyon](https://github.com/spencerlyon2)
* [Pablo Winant](http://www.mosphere.fr/)

The `gensys` and `csminwel` routines [`DSGE.gensys`](@ref) and
[`DSGE.csminwel`](@ref) are based on routines originally copyright
[Chris Sims](http://www.princeton.edu/~sims). The files are released here with
permission of Chris Sims under the BSD-3 [License](@ref).

The `kalman_filter` routine is loosely based on a version of the Kalman filter
algorithm originally copyright Federal Reserve Bank of Atlanta and written by
[Iskander Karibzhanov](http://karibzhanov.com). The files are released here with
permission of the Federal Reserve Bank of Atlanta under the BSD-3
[License](@ref).
