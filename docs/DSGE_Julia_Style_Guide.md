# DSGE Julia Style Guide

## Intro 

This document lists Julia coding recommendations consistent with best
practices in the software development community. The recommendations are
based on guidelines for other languages collected from a number of
sources and on personal experience. These guidelines are written with
the FRBNY DSGE code in mind. All pull requests submitted should follow these general style guidelines.


## Naming conventions

Emphasize readability! Our goal is for the code to mimic
the mathematical notation used in FRBNY DSGE papers as closely as possible.

### Variables

1. The names of variables should document their meaning or
use. Variables with a large scope should have especially meaningful
names. Variables with a small scope can have short names.

2. Exhibit consistency with the existing codebase.

3. Variables with economic or statistical significance should take
unicode syntax or imitate LaTeX syntax.  For example, ρ often
signifies the coefficient on an AR(1) process, and σ usually stands
for standard deviation. Parameters in the text should keep the same
symbol in the code (e.g. α in the code is the same α as in [this
paper](http://www.newyorkfed.org/research/staff_reports/sr647.html),
and takes on it usual significance as the capital share in a
Cobb-Douglas output function.


4. *Underscores:* 

  -Variable names should be in lower case, using underscores to
separate parts of a compound variable name. For example,
`steady_state` and `equilibrium_conditions` are two fields in the
`Model990()` object that follow this convention. Also notice that,
though the words could be shortened, they are spelled out for maximum
clarity.

  -Consistent with (3), underscores can and should be used when the
variable refers to a mathematical object that has a subscript. (In
this case, we are imitating LaTeX syntax.) For example, $r_m$ in LaTeX
should be represented by the variable `r_m`.

  -If the mathematical object has multiple subscripts, for example
$x_{i,j}$, simply concatenate the subscripts into `x_ij`.

  -If the object has superscripts as well as subscripts, for example $y^f_t$,
separate the superscripts with an underscore and place them first:
`y_f_t`.

  -For compatibility with existing code, variables with numeric
subscripts should exclude the underscore (e.g. `G0`, `ψ1`).


5. *Suffixes:*

  -Time: Consistent with the previous bullet points, the suffix `_t` as in
 `x_t` signifies the value of `x` at time `t`. The suffix `_t1`
 signifies the value of `x` at time `t-1`.

  -Shocks: The suffix `_sh` refers to a model shock.

6. *Prefixes:* 

  -The prefix `eq_` refers to an equilibrium condition.
 
  -The prefix `E` refers to an expecational shock.
 
  -The prefix `num_` should be used for variables representing the
number of objects (e.g. `num_parameters` or `num_anticipated_shocks`)
use the suffix `s` as is natural in spoken language.

  - Observables with the prefix `g` refer to growth rates.

7. Negative Boolean variable names should be avoided. A problem arises
when such a name is used in conjunction with the logical negation
operator as this results in a double negative. It is not immediately
apparent what `!isNotFound` means.  Use `isFound`. Avoid `isNotFound`.

8. Matrices that have mathematical significance (e.g. inputs to
Gensys, or the matrices of the transition and measurement equations)
should be upper case, as they are in mathematical notation (e.g. `TTT` or `YY`).




## Parameters


- A parameter vector is an object of user-defined type `parameter`
- A parameter vector of type `Parameters` collects more fundamental
  individual objects of type `Param`, which have fields

  - `Value`:   Float64
  - `IsFixed`: Logical
  - `Boundaries`
  - `PriorDistribution`
  - `TransformationType`: To go from model space to real line & vice versa

We define the following functions to act on objects of type parameter:

- `PriorDensity`: Look at value of parameter & its prior density; compute
- `Transform`:    To go from model space to real line
- `InvTransform`: To go from real line to model space

Parameter file for a particulare model looks like this, for all parameters:
```
# Define individual parameters as:
#
#   paraname = Param(Value, IsFixed, Bounds, PriorDistribution)
#
α = Param(0.2, false,[-1,1],     Beta(1,0.5))
β = Param(27,  true, [-Inf,Inf], Gamma(1,0.5))

# Collect parameters in vector
θ = Parameters(α, β, ...)
```

## Defining Indices (deprecated)

We have several functions for defining model indices, all collected into a
"ModelInds990" file.

All functions take a name like "π_t" or "rm" and give back an index.
Functions are:

- `endo`: For endogenous states
- `exo`:  Exogenous shocks
- `exp`:  Expectation shocks
- `eq`:   Equation indices
- `obs`:  Indices of named observables to use in measurement equation

They all look like this:
```julia
  function exo(name)
    names = ["π_sh", "rm_sh", "jerry", "george", "elaine", "kramer"]
    return find(map(nm -> (nm == name), names))
  end
```
- Since we don't care about the number, we only have to define the names.
- In this setup, adding states is easier, because we don't have to
  increment the index numbers of _everything_ when we add states.
- Super-automatic and less error prone; code focuses on the names just
  like we do.

## Equilibrium Conditions

A model-specific function of parameters and indices. Should look very
similar to our current code. Example
```julia
function eqcond990(θ::Parameters, endo::EndoStates, exo::ExoShocks, exp::ExpShocks, eq::Equations)

  Γ0[eq["mp"], endo["R_t"]) = 1;
  Γ1[eq["mp"], endo["R_t"]) = θ.ρ;
  Γ0[eq["mp"], endo["π_t"]) = -θ.Ψ_1;
  etc.

  return Γ0, Γ1, Ψ, Π

end
```
Measurement Equation will be very similar, taking parameters, model
indices, and data.

## Further organize model logic into types

-	Create `State` type, with fields 
    -	`name`
    -	`index`
    -	`description`
    -	`adjustmenttype` (e.g. cum_for in matlab code)
    -	`constant` (e.g. C_ss in matlab code)
-	`States` abstract type
-	`EndoStates <: States` abstract type
-	`EndoStates990 <: EndoStates` concrete type
-	`ExoShocks <: States` abstract type
-	`ExoShocks990 <: ExoShocks` concrete type
-	`ExpShocks <: States` abstract type
-	`ExpShocks990 <: ExpShocks` concrete type
-	`Equation` type
-	`Equations` abstract type
-	`Equations990 <: Equations` concrete type

Add additional functionality to the model types
- We should be able to iterate through `Param`s in a `Parameters` type:
```julia
assert(isa(θ, Parameters))
for α in θ
    println("$α")
end
```
- We should potentially be able to iterate through `State`s in instances of type `States`.
- We should be able to use the get functionality in instances of type `States` or `Equations`:
```julia
julia> eq["euler"]
1

julia> endo["c_t"]
1

julia> Γ0[eq["euler"], endo["c_t"]] = 1
```

## Defining a Model

A model is a user defined type with the following fields

- `Parameters`
- `EndoStates`
- `ExoShocks`
- `ExpShocks`
- `Equations`
- `Eqcond` (function)
- `Measurement` (function)

That's essentially enough to define the entire model structure.

We then build functions to act on a Model Type

- `dsgesolv`: Take a model type, use its parameters, indices, and eqconds to run through gensys and get TTT, RRR, CCC
- `Estimate`: Call dsgesolv, run gibb, etc.


## Building From There

Everything else will look like our current code, and will be more a less
a direct port of functions (Gensys, csminwel, kalman filter, etc.).

The main difference will be that rather than pass "TTT", "RRR", "YY" as
disjoint variables, we pass an object of "Model type" that encodes the
stucture, and _unpack_ that structure, either within functions or when
we pass arguments to a function.

We can also easily query key information about the model (whether
parametrs are fixed, what the prior distributions are, the index of
state "π_t"), since _one, single object_ of type "Model" will contain
all the relevant information.

