# DSGE.jl

The *DSGE.jl* package implements the FRBNY DSGE model and provides
general code to estimate many user-specified DSGE models. The package is introduced in the Liberty Street Economics blog post [The FRBNY
DSGE Model Meets Julia](http://libertystreeteconomics.newyorkfed.org/2015/12/the-frbny-dsge-model-meets-julia.html).

This Julia-language implementation mirrors the MATLAB code
included in the Liberty Street Economics blog post
[The FRBNY DSGE Model Forecast](http://libertystreeteconomics.newyorkfed.org/2015/05/the-frbny-dsge-model-forecast-april-2015.html).

FRBNY is currently working on extending the code to include forecasts and other features. Extensions of the DSGE model code may be released in the future at the discretion of FRBNY.

# Table of contents

```@contents
Depth = 4
```

# Acknowledgements
Developers of this package at [FRBNY](https://www.newyorkfed.org/research) include

* [Pearl Li](https://github.com/pearlzli)
* [Erica Moszkowski](https://github.com/emoszkowski)
* [Micah Smith](https://github.com/micahjsmith)

Contributors to this package at [QuantEcon](http://quantecon.org) include

* [Zac Cranko](https://github.com/ZacCranko)
* [Spencer Lyon](https://github.com/spencerlyon2)
* [Pablo Winant](http://www.mosphere.fr/)

The `gensys` and `csminwel` routines in [gensys.jl](src/solve/gensys.jl) and
[csminwel.jl](src/estimate/csminwel.jl) are based on routines originally
copyright [Chris Sims](http://www.princeton.edu/~sims). The files are released
here with permission of Chris Sims under the BSD-3 [license](LICENSE).

The `kalman_filter` routine is loosely based on a version of the
Kalman filter algorithm originally copyright Federal Reserve Bank of Atlanta
and written by [Iskander Karibzhanov](http://karibzhanov.com). The files are
released here with permission of the Federal Reserve Bank of Atlanta under the
BSD-3 [license](LICENSE).
