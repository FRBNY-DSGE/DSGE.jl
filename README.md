# New York Fed DSGE Model (Version 1002)
![Build Status](https://github.com/FRBNY-DSGE/DSGE.jl/workflows/build/badge.svg?branch=main)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://frbny-dsge.github.io/DSGE.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://frbny-dsge.github.io/DSGE.jl/latest)
[![codecov](https://codecov.io/gh/FRBNY-DSGE/DSGE.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/FRBNY-DSGE/DSGE.jl)

The `DSGE.jl` package implements the New York Fed dynamic stochastic general equilibrium (DSGE) model and provides general code to estimate many user-specified DSGE models. The package is introduced in the Liberty Street Economics blog post
[The FRBNY DSGE Model Meets Julia](http://libertystreeteconomics.newyorkfed.org/2015/12/the-frbny-dsge-model-meets-julia.html).
(We previously referred to our model as the "FRBNY DSGE Model.")

This Julia-language implementation mirrors the MATLAB code included in the
Liberty Street Economics blog post
[The FRBNY DSGE Model Forecast](http://libertystreeteconomics.newyorkfed.org/2015/05/the-frbny-dsge-model-forecast-april-2015.html).

Documentation for the *code* can be accessed by clicking on the `docs|dev` button above. For documentation about the most recent *model version*, read this [pdf](https://github.com/FRBNY-DSGE/DSGE.jl/blob/main/docs/DSGE_Model_Documentation_1002.pdf).

The New York Fed DSGE team is currently extending the code to solve and estimate heterogeneous agent models. Filtering and smoothing algorithms are available in the registered package [StateSpaceRoutines.jl](https://github.com/FRBNY-DSGE/StateSpaceRoutines.jl).
An implementation of Sequential Monte Carlo (SMC) sampling, used for the estimation of DSGE models, can be found in the registered package [SMC.jl](https://github.com/FRBNY-DSGE/SMC.jl). The foundational `AbstractModel` type, from which the `AbstractDSGEModel` type derives, is defined in the registered package [ModelConstructors.jl](https://github.com/FRBNY-DSGE/ModelConstructors.jl).

Further extensions of the DSGE model code may be released at the discretion of the New York Fed.

## Installation

`DSGE.jl` is a registered Julia package in the [`General`](https://github.com/JuliaRegistries/General) registry. To install it, open your Julia REPL, type `]` to enter the package manager, and run

```julia
pkg> add DSGE
```

If you use any code that loads data (e.g. the example script `run_default.jl` and `make_packet.jl`), then you need make sure you have a FRED API key by following these [instructions for the FredData.jl package](https://github.com/micahjsmith/FredData.jl).

If you are using Windows OS and you encounter the error `AssertionError: length(dirs) == 1`, please see this [issue](https://github.com/JuliaLang/Pkg.jl/issues/1943). Additionally, please do not run the `plot.jl` test if you are using Windows OS
because the generated output will violate the default filename length restriction on Windows. If you want to run this test, then
you need to enable [long paths](https://docs.microsoft.com/en-us/windows/win32/fileio/naming-a-file#enable-long-paths-in-windows-10-version-1607-and-later).

*Note we do not test our code in Windows OS, so we cannot guarantee the code works properly in Windows.*

## Versioning

`DSGE.jl` is currently compatible with Julia `v1.x` (as of `v1.1.6`).

To use `DSGE.jl` with Julia `v0.7`, please check out tag `0.8.1`. To do this, click on the drop-down menu that reads `branch:main` on the left-hand side of the page. Select `tags`, then `v0.8.1`.  If you've already cloned the repo, you can simply run `git checkout v0.8.1`.

To use `DSGE.jl` with Julia `v0.6`, please check out tag `0.4.1`.

## Precompilation

The `DSGE.jl` package is not precompiled by default because when running code in parallel, we want to re-compile
the copy of `DSGE.jl` on each processor to guarantee the right version of the code is being used. If users do not
anticipate using parallelism, then users ought to change the first line of `src/DSGE.jl` from

```
isdefined(Base, :__precompile__) && __precompile__(false)
```

to

```
isdefined(Base, :__precompile__) && __precompile__(true)
```

## Citing DSGE.jl

DSGE.jl (Version 1.2.1)[Source code], https://github.com/PSLmodels/DSGE.jl

Disclaimer
------
Copyright Federal Reserve Bank of New York. You may reproduce, use, modify, make derivative works of, and distribute and this code in whole or in part so long as you keep this notice in the documentation associated with any distributed works. Neither the name of the Federal Reserve Bank of New York (FRBNY) nor the names of any of the authors may be used to endorse or promote works derived from this code without prior written permission. Portions of the code attributed to third parties are subject to applicable third party licenses and rights. By your use of this code you accept this license and any applicable third party license.

THIS CODE IS PROVIDED ON AN "AS IS" BASIS, WITHOUT ANY WARRANTIES OR CONDITIONS OF ANY KIND, EITHER EXPRESS OR IMPLIED, INCLUDING WITHOUT LIMITATION ANY WARRANTIES OR CONDITIONS OF TITLE, NON-INFRINGEMENT, MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, EXCEPT TO THE EXTENT THAT THESE DISCLAIMERS ARE HELD TO BE LEGALLY INVALID. FRBNY IS NOT, UNDER ANY CIRCUMSTANCES, LIABLE TO YOU FOR DAMAGES OF ANY KIND ARISING OUT OF OR IN CONNECTION WITH USE OF OR INABILITY TO USE THE CODE, INCLUDING, BUT NOT LIMITED TO DIRECT, INDIRECT, INCIDENTAL, CONSEQUENTIAL, PUNITIVE, SPECIAL OR EXEMPLARY DAMAGES, WHETHER BASED ON BREACH OF CONTRACT, BREACH OF WARRANTY, TORT OR OTHER LEGAL OR EQUITABLE THEORY, EVEN IF FRBNY HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGES OR LOSS AND REGARDLESS OF WHETHER SUCH DAMAGES OR LOSS IS FORESEEABLE.
