# Notes for DSGE.jl Contributors

## DSGE Julia Style Guide

### Intro 

This document lists Julia coding recommendations consistent with best
practices in the software development community. The recommendations are
based on guidelines for other languages collected from a number of
sources and on personal experience. These guidelines are written with
the FRBNY DSGE code in mind. All pull requests submitted should follow these
general style guidelines.

### Naming conventions

Emphasize readability! Our goal is for the code to mimic the mathematical
notation used in FRBNY DSGE papers as closely as possible.

1. The names of variables should document their meaning or
use. Variables with a large scope should have especially meaningful
names. Variables with a small scope can have short names.

2. Exhibit consistency with the existing codebase.

3. Modules and type names in UpperCamelCase

4. Variable names in snake_case

Variable names should be in lower case, using underscores to
separate parts of a compound variable name. For example,
`steady_state` and `equilibrium_conditions` are two fields in the
`Model990()` object that follow this convention. Also notice that,
though the words could be shortened, they are spelled out for maximum
clarity.

5. The prefix `num_` should be used for variables representing the
number of objects (e.g. `num_parameters` or `num_anticipated_shocks`)
use the suffix `s` as is natural in spoken language.

6. Negative Boolean variable names should be avoided. A problem arises
when such a name is used in conjunction with the logical negation
operator as this results in a double negative. It is not immediately
apparent what `!isNotFound` means.  Use `isFound`. Avoid `isNotFound`.

7. Named constants can be all uppercase using underscore to separate words:
`MAX_ITERATIONS`

8. Naming mathematical objects

Variables with mathematical significance should use 
unicode characters and imitate LaTeX syntax.  For example, ρ should be used to  
name the autocorrelation coefficient in an AR(1) process, and σ should be used
to name standard deviation. Parameters in the text should keep the same
symbol in the code (e.g. α in the code is the same α as in [this
paper](http://www.newyorkfed.org/research/staff_reports/sr647.html),
and takes on it usual significance as the capital share in a
Cobb-Douglas output function.

- General conventions
  - Underscores can and should be used when the variable refers to a
    mathematical object that has a subscript. (In this case, we are
    imitating LaTeX syntax.) For example, $r_m$ in LaTeX should be
    represented by the variable `r_m`.
  - If the mathematical object has multiple subscripts, for example $x_{i,j}$,
    simply concatenate the subscripts: `x_ij`.
  - If the object has superscripts as well as subscripts, for example
    $y^f_t$, separate the superscripts with an underscore and place them
    first: `y_f_t`.
  - For compatibility with existing code, variables with numeric subscripts
    should exclude the underscore: `G0`, `ψ1`.
  - Matrices that have mathematical significance (e.g. the matrices of the
    transition and measurement equations) should be upper case, as they are
    in mathematical notation, and can repeat the letter to avoid collisions:
    `TTT` or `YY`.
  - Symbols such as overbars (which indicate mean values) and  tildes (which indicate
    log-deviations from the steady state) are written using a 3- or
    4-letter abbreviation immediately after the variable they modify:
    `kbar_t`, `ztil` (z tilde).
  - Stars indicating steady-state variables are included as
    subscripts: `π_star_t`
- Suffixes
  - Time: Consistent with the previous bullet points, the suffix `_t` as in
    `x_t` signifies the value of `x` at time `t`. The suffix `_t1`
    signifies the value of `x` at time `t-1`.
  - Shocks: The suffix `_sh` refers to a model shock.
- Prefixes
  - The prefix `eq_` refers to an equilibrium condition.
  - The prefix `E` refers to the expectation operator.
  - Observables with the prefix `g` refer to growth rates.

### Code Formatting Guidelines

- Indent 4 spaces
- Wrap lines at 92 characters
- Use whitespace to enhance readability
- No trailing whitespace
