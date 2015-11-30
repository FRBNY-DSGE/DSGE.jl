# Technical Appendix
*Zac Cranko, Pearl Li, Spencer Lyon, Erica Moszkowski, Micah Smith, Pablo Winant*

The FRBNY DSGE model is a relatively large New Keynesian model augmented with
financial frictions and a variety of innovations. Here at the Fed, we use it
both for forecasting and policy analysis. Research using this model includes
looking at the dynamics of inflation during the great recession, the effects of
forward guidance, and more. When we were approached by the folks at QuantEcon
about a possible collaboration, we jumped at the idea, as it would give us an
opportunity to rework our code in an arguably faster language, redesign it from
the ground up, and release it open source for the benefit of the community.
Julia was the language of choice, recommended by the QuantEcon group for its
high performance and suitability for this breed of technical computing. [Any
additional commentary about our collaboration from the QuantEcon folks would be
great!] In this post, we’ll discuss our experiences redesigning our code from
the ground up, the resulting performance changes, and the challenges we faced
working with such a young language.  We created a suite to assess the
performance of our Julia code, relative to our MATLAB code. We focus on both
the core functions used in solving and estimating the model, as well as on
longer-running routines of greater scope. These tests were conducted on a
computing cluster running RHEL 6.7 on Intel Xeon E5-2697 v2 (2.70Ghz) cores.

Table 1: Benchmark times relative to MATLAB (smaller is better)

| Test                 | MATLAB (14a) | Julia (v0.4.0) |
| -------------------- | :----------: | :------------: |
| `gensys`             | 1.00         | 0.17           |
| `solve`              | 1.00         | 0.09           |
| `kalman_filter`      | 1.00         | 0.75           |
| `posterior`          | 1.00         | 0.26           |
| `csminwel`           |              |                |
| `hessian`            | 1.00         | 0.23           |
| `metropolis_hastings`| 1.00         | 0.11           |

We ultimately achieve an increase of speed that reduces running time to 1/10th
to 3/4th that of the MATLAB code. This increase is very significant for our
purposes, as we reduce the running time of the "full shebang" (optimization,
hessian calculation, and Metropolis-Hastings sampling) by a substantial amount.
How much of this increase is due to native performance adventures of Julia, and
how much is simply due to the improvements in design that came from rebuilding
this project from the ground up? It is of course difficult to say. Our MATLAB
code suffers from inefficiencies due to its long, cumulative development; these
have been addressed in our Julia package. If we are most interested in
isolating differences to the languages themselves, our code to compute the
model solution with gensys and apply the Kalman filter has relatively little
redesign and optimization as compared to the MATLAB code and provides the most
comparable, though still imperfect, measurements of performance. This reduction
of 1/5th to 3/4th could therefore be taken as a reasonable estimate as to
Julia's advantage in this one arena of computation.

## Code Improvements

This section discusses how we used some of Julia’s features to improve upon our
old code. Julia’s robust and flexible type system has allowed us to improve our
code’s performance and clarity in several fundamental ways. First, Julia
supports an object-oriented programming paradigm that lends itself naturally to
our DSGE model.  *DSGE.jl* centers around a model object, which stores all
information associated with the model – including a fairly large number of
parameters, priors, states, equilibrium conditions, computational settings, and
flags – in one place.  By simply passing the model object as an argument to any
function, the function has access to all of the model’s fields. By comparison,
our old MATLAB code stored all variables directly the workspace – an approach
that scaled poorly as model specifications become more and more complex. To
illustrate just how unwieldy our MATLAB code was, many of our function calls
required more than 20 arguments (a nightmare for usage and human-readability).

Furthermore, an object-oriented approach allows us to take advantage of method
dispatch in Julia by defining different model types for different model
specifications. As detailed in *DSGE.jl*’s README file [LINK], changes to the
model’s equilibrium conditions and measurement equation are referred to as
changes in a model "specification."  In the Julia code, model specifications
have a 1:1 correspondence with types: every specification is implemented by a
Julia type. Where necessary, the same function can be defined differently for
different model types (for example, the augment_states function adds additional
states to the model’s transition matrices after it is solved. We can pass any
model object to augment_states, and Julia ensures that the proper version of
the function will be executed). Because MATLAB does not have multiple dispatch,
any functions that required different behavior for different model
specifications included many switch statements and if clauses – another
approach that made our MATLAB code cumbersome and clunky. Of course, multiple
dispatch improves our code’s performance more generally: type declarations
allow many functions to be precompiled, reducing the need for type inference at
runtime.

It is easy to see that all model types constructed for use with DSGE.jl will be
closely related: they will have the same fields, and are passed to the same
methods.  If it sounds to you like we have an implicit interface here, you’d be
right. Rather than implementing each object as a standalone type, we define an
abstract type, `AbstractDSGEModel`, to serve as the parent for all model types.
Because most of our routines are not model-specific, we need only define them
once (with an argument of type AbstractDSGEModel) and Julia’s type promotion
takes care of the rest. These functions expect the model object to have certain
fields, and for those fields to have certain types. In this way, although Julia
does not have a mechanism for explicit interface definitions (as Java does, for
example), its type system enforces an implicit DSGE model interface. With a
clear interface in place, running new model specifications using DSGE.jl is
relatively straightforward. (See here [LINK to “Editing or Extending a Model”
section of README] for detailed instructions).  We also define model parameters
[LINK TO README] as subtypes of a common abstract type AbstractParameter. This
allows us to abstract to one notion of a model parameter, while implementing
different kinds of parameters in different ways (for example, steady-state
parameters are calculated as a function of other parameters, so have different
required fields).

We use parameterized types to increase the flexibility of our codebase in a
number of ways. For example, we are now able to think of all computational
settings, such as flags, data vintages, and the number of iterations for which
we run the Metropolis-Hastings sampling algorithm, as objects of type Setting.
We can have settings that are Booleans, Integers, Floats, and Strings, and
store them centrally in a single, type-stable dictionary within the model
object.

Julia's JIT compilation also provides significant performance boosts in some
areas. For example, we allow in our model a variable number of anticipated
monetary policy shocks, beginning in 2007Q1, that we use to treat the zero
lower bound. In our MATLAB code, we suffer some dynamic code generation to
implement this feature. Julia's compile-time evaluation of such statements
eliminates this performance hit. Granted, there may be better solutions to our
problem in both languages, but similar situations involving code generation are
easily addressed in Julia. Then, of course, the presence of compilation itself
(as opposed to purely interpreted code) provides further performance
improvement, of which the Julia community has documented widely.

Finally, we have found that a number of Julia features make working with
DSGE.jl much more pleasant than working with our old codebase. Julia’s testing
infrastructure has significantly improved our development workflow. Unicode
support means that code can correspond more closely to actual model equations,
reducing the headache associated with translating between math-speak and
code-speak. Operator overloading and user-defined syntax makes it easy to
access different values (for example, we can access the value of model m’s
parameter α using the simple syntax m[:α], and we use the syntax m <= Foo  to
add settings and parameters to the model object). Inline markdown documentation
also helps improve the developer experience.

Possible improvements to DSGE.jl are many and varied. We may consider using
alternative optimization routines to improve speed. Ultimately, powerful
metaprogramming support (such as a macro to parse equilibrium conditions) would
allow user to specify model equations almost literally. We welcome such
improvements to the existing code from the community – see [LINK
CONTRIBUTING.md].

## Challenges
Converting the FRBNY DSGE model from MATLAB, a mature and well-supported
language, to an extremely young language like Julia involves no shortage of
challenges. Significant changes to the Julia language itself are introduced in
rapid succession, and using DSGE.jl with a new Julia version inevitably floods
the user’s screen with deprecation warnings. There is significant difficulty in
finding written resources on the language in addition to the Julia Manual.
Google searches frequently return discussions in Issues, which are unhelpful to
elementary users and can be actively misleading at times. 

Differences between the behavior of MATLAB and Julia’s core linear algebra
libraries lead to many roadblocks in the development of DSGE.jl. Julia uses
BLAS functions, which use multithreading, for some linear algebra functions.
Using a different number of threads changes the results of matrix decomposition
when the matrix is singular.   This indeterminacy caused significant problems
for our testing suite, both in comparing output matrices to MATLAB results and
in testing for reproducibility among Julia outputs. To retain as much testing
coverage as possible while still taking advantage of multithreading, we
implemented a separate test to test eigenvalue decomposition using a single
compute thread. [Maybe Pearl/Pablo can comment on ordered qz decomposition and
eigenvalue/singular value decomposition here?] We ran into most of these
problems while porting the model solution algorithm, gensys.  TFor example, the
MATLAB model code uses MATLAB’s default behavior of returning a possibly
complex QZ (generalized Schur) decomposition A = QSZ* and B = QTZ* with upper
triangular matrices S and T, while the default Julia behavior returns a real
decomposition with upper Hessenberg (with blocks on the diagonals) S and T. 

Finally, dealing with a recently introduced language can make it quite
difficult for new users to produce performant code.  The typical economist,
especially one coming from a MATLAB background, may be unfamiliar with the
nature and use of concepts like type stability, parametric types, and
preallocation. Julia’s profiler and debugger lack the flexibility of those in
MATLAB, and can make it difficult to identify the source of slow performance or
errors. However, we’re confident that Julia syntax should be familiar and
intuitive enough to be picked up easily by MATLAB programmers, who should be
able to jump right in to Julia without too much immediate overhead. [Zac, any
additions to this paragraph?]

## Disclaimer
This post reflects the experience of the authors with Julia and MATLAB and does
not represent an endorsement by the Federal Reserve Bank of New York or the
Federal Reserve System of any particular product or service. The views
expressed in this post are those of the authors and do not necessarily reflect
the position of the Federal Reserve Bank of New York or the Federal Reserve
System. Any errors or omissions are the responsibility of the authors.
