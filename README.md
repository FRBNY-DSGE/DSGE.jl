# New York Fed DSGE Model (Version 1002)
[![Build Status](https://travis-ci.org/FRBNY-DSGE/DSGE.jl.svg)](https://travis-ci.org/FRBNY-DSGE/DSGE.jl)
[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://frbny-dsge.github.io/DSGE.jl/latest)

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

The New York Fed DSGE team is currently extending the code to solve and estimate
heterogeneous agent models.
An implementation of Sequential Monte Carlo (SMC) sampling can be found in the registered package [SMC](https://github.com/FRBNY-DSGE/SMC.jl).
Further extensions of the DSGE model code may be released in the future at the discretion
of the New York Fed.

## Installation

`DSGE.jl` is a registered Julia package. To install it, open your Julia REPL, type `]` (enter package manager), and run

```julia
pkg> add DSGE
```
## Versioning

`DSGE.jl` is currently compatible with 0.7. Official compatibility with 1.0 is imminent.

To use `DSGE.jl` with Julia version 0.6, please check out tag
0.4.1. To do this, click on the drop-down menu that reads `branch:
master` on the left-hand side of the page. Select `tags`, then
`v0.4.1`.  If you've already cloned the repo, you can simply run
`git checkout v0.4.1`.

Disclaimer
------
Copyright Federal Reserve Bank of New York. You may reproduce, use, modify, make derivative works of, and distribute and this code in whole or in part so long as you keep this notice in the documentation associated with any distributed works. Neither the name of the Federal Reserve Bank of New York (FRBNY) nor the names of any of the authors may be used to endorse or promote works derived from this code without prior written permission. Portions of the code attributed to third parties are subject to applicable third party licenses and rights. By your use of this code you accept this license and any applicable third party license.

THIS CODE IS PROVIDED ON AN "AS IS" BASIS, WITHOUT ANY WARRANTIES OR CONDITIONS OF ANY KIND, EITHER EXPRESS OR IMPLIED, INCLUDING WITHOUT LIMITATION ANY WARRANTIES OR CONDITIONS OF TITLE, NON-INFRINGEMENT, MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, EXCEPT TO THE EXTENT THAT THESE DISCLAIMERS ARE HELD TO BE LEGALLY INVALID. FRBNY IS NOT, UNDER ANY CIRCUMSTANCES, LIABLE TO YOU FOR DAMAGES OF ANY KIND ARISING OUT OF OR IN CONNECTION WITH USE OF OR INABILITY TO USE THE CODE, INCLUDING, BUT NOT LIMITED TO DIRECT, INDIRECT, INCIDENTAL, CONSEQUENTIAL, PUNITIVE, SPECIAL OR EXEMPLARY DAMAGES, WHETHER BASED ON BREACH OF CONTRACT, BREACH OF WARRANTY, TORT OR OTHER LEGAL OR EQUITABLE THEORY, EVEN IF FRBNY HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGES OR LOSS AND REGARDLESS OF WHETHER SUCH DAMAGES OR LOSS IS FORESEEABLE.
