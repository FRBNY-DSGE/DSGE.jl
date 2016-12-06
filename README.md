# FRBNY DSGE Model (Version 990.2), Sequential Monte Carlo Development Branch
[![Build Status](https://travis-ci.org/FRBNY-DSGE/DSGE.jl.svg)](https://travis-ci.org/FRBNY-DSGE/DSGE.jl)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://FRBNY-DSGE.github.io/DSGE.jl/stable)
[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://FRBNY-DSGE.github.io/DSGE.jl/latest)

The *DSGE.jl* package implements the FRBNY DSGE model and provides
general code to estimate many user-specified DSGE models. The package
is introduced in the Liberty Street Economics blog post [The FRBNY DSGE Model Meets Julia](http://libertystreeteconomics.newyorkfed.org/2015/12/the-frbny-dsge-model-meets-julia.html).

This Julia-language implementation mirrors the MATLAB code included in
the Liberty Street Economics blog post [The FRBNY DSGE Model Forecast](http://libertystreeteconomics.newyorkfed.org/2015/05/the-frbny-dsge-model-forecast-april-2015.html).

This branch of the FRBNY *DSGE.jl* package also implements Sequential
Monte Carlo (SMC) sampling as an alternative to Metropolis Hastings
Markov Chain Monte Carlo sampling. This additional code is based upon
Edward Herbst and Frank Schorfheide's 2014 paper "Sequential Monte
Carlo Sampling for DSGE Models" and the code accompanying their book
*Bayesian Estimation of DSGE Models*. Currently, FRBNY's
implementation of SMC works on the small-scale New Keynesian DSGE
model presented in Sungbae An and Frank Schorfheide's 2007 paper
"Bayesian Analysis of DSGE Models". FRBNY is currently working on
extending the code so that SMC may be used with medium-scale DSGE
models. This and other extensions of the DSGE model code may be
released in the future at the discretion of FRBNY.

## Usage

Once you've followed the installation instructions for the main
`DSGE.jl` package, you can switch the the SMC branch by opening your
terminal, navigating to the top level `DSGE.jl` directory, and running
`git checkout smc`. To return to the original *DSGE.jl* package, run
`git checkout master`. The file `test_smc.jl` in
`DSGE\docs\examples\` is the appropriate place for a new user to
start. It initializes a new model, sets some SMC-related parameters,
runs SMC, and generates a Latex document documenting parameter moments. 

Currently, SMC works only on the small-scale An Schorfheide model. On
medium-scale models such as model 990, SMC frequently fails. Current
work on this branch is focused on addressing this issue, as well as
implementing blocked parameter sampling and an adaptive tuning
schedule (as discussed in "Sequential Monte Carlo Sampling for DSGE
Models"). Comments and suggestions are welcome.