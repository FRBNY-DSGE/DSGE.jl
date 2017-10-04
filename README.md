# New York Fed DSGE Model (Version 1002)

[![Build Status](https://travis-ci.org/FRBNY-DSGE/DSGE.jl.svg)](https://travis-ci.org/FRBNY-DSGE/DSGE.jl)
[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://FRBNY-DSGE.github.io/DSGE.jl/latest)

The *DSGE.jl* package implements the New York Fed DSGE model and provides
general code to estimate many user-specified DSGE models. The package is
introduced in the Liberty Street Economics blog post
[The FRBNY DSGE Model Meets Julia](http://libertystreeteconomics.newyorkfed.org/2015/12/the-frbny-dsge-model-meets-julia.html).
(We previously referred to our model as the "FRBNY DSGE Model".)

This Julia-language implementation mirrors the MATLAB code included in the
Liberty Street Economics blog post
[The FRBNY DSGE Model Forecast](http://libertystreeteconomics.newyorkfed.org/2015/05/the-frbny-dsge-model-forecast-april-2015.html).

For the latest documentation on the *code*, click on the docs|latest button
above. Documentation for the most recent *model version* is available
[here](https://github.com/FRBNY-DSGE/DSGE.jl/blob/master/docs/DSGE_Model_Documentation_1002.pdf).

The New York Fed DSGE team is currently working on extending the code to include
forecasts and other features. Ongoing work on implementing Sequential Monte
Carlo (SMC) sampling can be found on the `smc` branch. Further extensions of the
DSGE model code may be released in the future at the discretion of the New York
Fed.

## Installation

`DSGE.jl` is a registered Julia package. To install it, open your Julia REPL and run

```julia
julia> Pkg.add("DSGE")
```
# Versioning

To use `DSGE.jl` with Julia version 0.4, please check out tag
0.2.0. To do this, click on the drop-down menu that reads `branch:
master` on the left-hand side of the page. Select `tags`, then
`v0.2.0`.  If you've already cloned the repo, you can simply run
`git checkout v0.2.0`.
