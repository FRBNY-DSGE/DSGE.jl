# Herbst-Schorfheide Sequential Monte Carlo Implementation and Replication
[![Build Status](https://travis-ci.org/FRBNY-DSGE/DSGE.jl.svg)](https://travis-ci.org/FRBNY-DSGE/DSGE.jl)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://FRBNY-DSGE.github.io/DSGE.jl/stable)
[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://FRBNY-DSGE.github.io/DSGE.jl/latest)


This branch of the FRBNY *DSGE.jl* package implements Sequential Monte
Carlo (SMC) sampling as an alternative to Metropolis Hastings Markov
Chain Monte Carlo sampling. The SMC algorithm implemented here is
based upon Edward Herbst and Frank Schorfheide's paper ["Sequential
Monte Carlo Sampling for DSGE
Models"](http://dx.doi.org/10.1002/jae.2397) and the code accompanying
their book *Bayesian Estimation of DSGE Models*. More information and
the original Matlab scripts that this code replicates can be found at
Frank Schorfheide's
[website](https://sites.sas.upenn.edu/schorf/pages/bayesian-estimation-dsge-models).
Currently, FRBNY's implementation of SMC works on the small-scale New
Keynesian DSGE model presented in Sungbae An and Frank Schorfheide's
paper ["Bayesian Analysis of DSGE Models"](
http://dx.doi.org/10.1080/07474930701220071).  FRBNY is currently
working on extending the code so that SMC may be used with
medium-scale DSGE models. This and other extensions of the DSGE model
code may be released in the future at the discretion of
FRBNY. Comments and suggestions are welcome, and best submitted as
either an issue or a pull request to this branch.

## Background

The *DSGE.jl* package implements the FRBNY DSGE model and provides
general code to estimate many user-specified DSGE models. The package
is introduced in the Liberty Street Economics blog post [The FRBNY
DSGE Model Meets
Julia](http://libertystreeteconomics.newyorkfed.org/2015/12/the-frbny-dsge-model-meets-julia.html).

This Julia-language implementation mirrors the MATLAB code included in
the Liberty Street Economics blog post [The FRBNY DSGE Model
Forecast](http://libertystreeteconomics.newyorkfed.org/2015/05/the-frbny-dsge-model-forecast-april-2015.html).

## Usage

Once you've followed the installation instructions for the main
`DSGE.jl` package, you should you also have Git installed on your
machine. If you do not, you can download it from
[https://git-scm.com/](https://git-scm.com/). Once Git is installed,
you can switch to the SMC branch by opening your
Julia REPL and running `Pkg.checkout("DSGE","smc")`. To return to the
original *DSGE.jl* package, run `Pkg.checkout("DSGE")` from
the REPL. The file `test_smc.jl` in `DSGE/docs/examples/` is the
appropriate place for a new user to start. It initializes a new model,
sets some SMC-related parameters, runs SMC, and generates a LaTeX
document documenting parameter moments. You can run this file by
running `include("$(Pkg.dir("DSGE"))/docs/examples/test_smc.jl")` from
the REPL. 
