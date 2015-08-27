# DSGE Julia Style Guide

## Intro 


The following objects define a model: 

- __Parameters__: Have values, bounds, fixed-or-not status, priors
- __States__: Collections of type `State` that map a name to an index.
  (e.g. "π_t" -> 1)
- __Equilibrium Conditions__: A function that takes parameters and model
  indices, then returns Γ0, Γ1, Ψ, and Π

These are enough to define the model structure. _Everything else_ is
essentially a function of these basics, and we can get to a forecast by
this chain:

- (Parameters + Model Indices + Eqcond Function) -> (TTT + RRR)
- (TTT + RRR + Data) -> Estimation
- (Estimation + TTT + RRR + Data) -> Forecast


## The Model Object

A big improvement of the Julia code over Matlab is that we can follow the object-oriented 
programming paradigm. Rather than having a lot of different variables that we associate with 
different model specifications floating around, we can define different model types and easily
create a new instance of any type (with, say, new data or different flags). 
Each variable is expressly tied to a particular instance of that type of model, so we can pass
the whole model object to any function we want and that function will have access to the particular 
state associated with that instance of the model.

A model is defined with the following fields, which are guaranteed by the AbstractDSGEModel



## Naming conventions

Emphasize readability! This is our opportunity to make the code as self-documenting and clear
as possible. Always air on the side of more descriptive, rather than less.

### Variables

1. The names of variables should document their meaning or use. 

2. Exhibit consistency with the existing codebase.

3. Variable names should be in lower case, using underscores to separate parts of a compound variable name. 
  - For example, `steady_state` and `parameters_fixed` are two fields in the `Model990()` object that follow 
this Julia convention. 
  - Acronyms, even if normally uppercase, should be lower case.

4. Variables with economic or mathematical significance should imitate LaTeX syntax. Underscores can and 
should be used when the variable refers to a mathematical object that has a subscript (as in LaTeX).
  - For example, $r_m$ should be represented by the variable `r_m`. If the mathematical object 
has multiple subscripts, for example x_{i,j}, simply concatenate the subscripts into `x_ij`. If the object has
superscripts as well, for example $y^f_t$, separate the superscripts with an underscore and place them first: 
`y_f_t`.
  - For compatibility with existing code, variables with numeric subscripts should exclude the underscore:
`G0`, `psi1`

5. Variables with statistical significance (punny!) should be unicode characters. 
  - For example, μ usually signifies a mean, and σ usually stands for standard deviation. 
  - However, model parameters are currently defined as symbols that reflect LaTeX syntax. For example, the 
symbol for the parameter ψ in the parameters vector is :psi.

6. Variables with a large scope should have meaningful names. Variables with a small scope can have short names.
  - In practice most variables should have meaningful names. The use of short names should be reserved for 
conditions where they clarify the structure of the statements. Scratch variables used for temporary storage or 
indices can be kept short. A programmer reading such variables should be able to assume that its value is not 
used outside a few lines of code. Common scratch variables for integers are `i`, `j`, `k`, `m`, `n` and for
doubles `x`, `y` and `z`

7. The prefix *num_* should be used for variables representing the number of objects.
  - Enhances the readability of the code
  - Use the suffix `s` as is natural in spoken language (e.g. `num_parameters`)

8. A convention on pluralization should be followed consistently.
  - A suggested practice is to make all variable names either singular or plural. 
  - *Having two variables with names differing only by a final letter *s* should be avoided.*  An acceptable 
  alternative for the plural is to use the suffix `Array` or `List` (eg `point`, `point_list`)

9. Consistent with programming convention, iterator variables should be named or prefixed with `i`, `j`, `k`, etc.
  ```julia
   for i = 1:num_parameters(m)
      ;
   end
  ```

10. Named constants can be all uppercase using underscore to separate words.
  – e.g. `MAX_ITERATIONS`

## The Param Type

- The model object contains an array of 
- fundamental individual objects of type `Param`, which have fields

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

We have five functions for defining model indices, all collected into a
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

