# The DSGE Matlab to Julia Transition: Improvements and Challenges
*Zac Cranko, Pearl Li, Spencer Lyon, Erica Moszkowski, Micah Smith, Pablo Winant*

*December 3, 2015*

The FRBNY DSGE model is a relatively large New Keynesian model augmented with
financial frictions and a variety of innovations. Here at the Fed, we use it
both for forecasting and policy analysis. Research using this model includes
looking at the dynamics of inflation during the great recession, the effects of
forward guidance, and more. When we were approached by the folks at QuantEcon
about a possible collaboration, we jumped at the idea, as it would give us an
opportunity to rework our code in an arguably faster language, redesign it from
the ground up, and release it open source for the benefit of the community.
Julia was the language of choice, recommended by the QuantEcon group for its
high performance and suitability for this breed of technical computing.

[Any
additional commentary about our collaboration from the QuantEcon folks would be
great!]

In this post, we’ll discuss our experiences redesigning our code from
the ground up, the resulting performance changes, and the challenges we faced
working with such a young language.  We created a suite to assess the
performance of our Julia code, relative to our MATLAB code. We focus on both
the core functions used in solving and estimating the model, as well as on
longer-running routines of greater scope. These tests were conducted on a
single core on an Intel® Xeon® E5-2697 v2 2.70GHz CPU running GNU/Linux:

Benchmark times relative to MATLAB (smaller is better)

| Test                 | MATLAB (14a) | Julia (0.4.0)  |
| -------------------- | :----------: | :------------: |
| `gensys`             | 1.00         | 0.17           |
| `solve`              | 1.00         | 0.09           |
| `kalman_filter`      | 1.00         | 0.75           |
| `posterior`          | 1.00         | 0.26           |
| `csminwel`           | 1.00         | 0.33           |
| `hessian`            | 1.00         | 0.23           |
| `metropolis_hastings`| 1.00         | 0.11           |

We ultimately achieve an increase of speed that reduces running time to 1/10th
to 3/4th that of the MATLAB code. This increase is very significant for our
purposes, as we substantially reduce the running time of the "full shebang"
(optimization, hessian calculation, and Metropolis-Hastings sampling).

How much of this increase is due to native performance adventures of Julia, and
how much is simply due to the improvements in design that came from rebuilding
this project from the ground up? It is of course difficult to say, and it is
important to emphasize that one cannot be sure what portion of the performance
increase can be attributed to inherent language features as opposed to design
differences. Indeed, our MATLAB code suffers from inefficiencies due to its
long, cumulative development, and support for a plethora of models and features.
These design issues have been largely addressed in our Julia package. If we are
most interested in isolating differences to the languages themselves, our code
to compute the model solution with gensys and apply the Kalman filter has
relatively little redesign and optimization as compared to the MATLAB code and
provides the most comparable, though still imperfect, measurements of
performance. This reduction of 1/5th to 3/4th could therefore be taken as a
first estimate of Julia's advantage in this single arena of computation.

## Code Improvements

Julia provides versatile language features that allow us to improve our code's
performance and clarity in several fundamental ways. First and foremost of these
is the highly integrated, robust, and flexible type system that lends itself
naturally to our DSGE model. At the center of the *DSGE.jl* package is the
*model object*. Here, one can store all information associated with the model –
including the numerous parameters, priors, states, equilibrium conditions,
computational settings, and flags – in one place.  By simply passing the model
object as an argument to any function, the function has access to all of the
model’s fields.  By comparison, our MATLAB code stored all variables
directly in the global workspace – an approach that scaled poorly as model
specifications become more and more complex. To illustrate just how unwieldy our
MATLAB code was, many of our function calls required more than 20 positional
arguments (a nightmare for usage and human-readability). Certainly, one could
approximate a "model object" in MATLAB by using its own object-oriented classes,
or by "bundling" model attributes into a `struct` or other data structure.
However, MATLAB classes are both relatively complicated and slower than
non-object implementations. And large `struct`s would be susceptible to
excessive dynamic field access.

Furthermore, a type-based approach allows us to take advantage of method
dispatch in Julia by defining different model types for different model
specifications. As detailed in the
[README file](https://github.com/FRBNY-DSGE/DSGE.jl/blob/master/README.md),
changes to the model’s equilibrium conditions and measurement equation are
referred to as changes in a model's "specification."  In the Julia code, model
specifications have a 1:1 correspondence with concrete types.  Where necessary,
a single *function* can have multiple *methods* defined, that are customized for
different model types.  (For example, the `augment_states` function adds
augments the model’s transition matrices after it has been solved.  We can
pass any model object to `augment_states`, and Julia ensures that the proper
version of the function will be executed.) [TODO Code to illustrate here?]
In MATLAB, any functions that required different behavior for different model
specifications included many switch statements and if clauses – another approach
that made our MATLAB code cumbersome and clunky. Of course, multiple dispatch
improves our code’s performance more generally: type declarations allow many
functions to be precompiled, reducing the need for type inference at runtime,
and significant code specialization is possible.

It is easy to see that all model types constructed for use with *DSGE.jl* are
closely related: they will have the same fields, and are passed to the same
methods.  If it sounds to you like we have an implicit interface here, you’re
right. Rather than implementing each object as a standalone type, we define an
abstract type, `AbstractModel`, to serve as the parent for all model types.
Because most of our routines are not model-specific, we need only define them
once (with an argument of type `AbstractModel`) and Julia’s dispatch system
takes care of the rest. These functions expect the model object to have certain
fields, and for those fields to have certain types. In this way, although Julia
does not have a mechanism for explicit interface definitions (as other
languages, such as Java), its type system enforces an implicit DSGE model
interface. With a clear interface in place, running new model specifications
using *DSGE.jl* is relatively straightforward. (See
[here](https://github.com/FRBNY-DSGE/DSGE.jl/blob/master/README.md#extending-or-editing-a-model)
for detailed instructions). We similarly define model parameters as subtypes of a
common abstract type `AbstractParameter`. This allows us to abstract to one notion
of a model parameter, while implementing different kinds of parameters in
different ways. Steady-state parameters, for example, are really functions of
other parameters, and have different implementation.

[TODO not really relevant to Improvement section?]
We use parameterized (generic) types to increase the flexibility of our codebase
in a number of ways. For example, we are now able to think of all computational
settings, such as flags, data vintages, and the number of samples to draw using
the Metropolis-Hastings algorithm, as instances of type `Setting`. We can have
settings of various data types and store them centrally in a single dictionary
within the model object.

Julia's JIT compilation also provides significant performance boosts in some
areas. For example, we allow a variable number of anticipated
monetary policy shocks, beginning in 2008Q4, that we use to treat the zero
lower bound. In our MATLAB code, we suffer some dynamic code generation to
implement this feature. Julia's compile-time evaluation of such statements
eliminates this performance hit. Granted, there may be better solutions to our
problem in both languages, but similar situations involving code generation are
easily addressed in Julia. Then, of course, the compilation itself (as opposed
to purely interpreted code) provides further performance improvement.

We have found that a number of Julia features make working with
*DSGE.jl* simply more pleasant and user-friendly than working with our old
codebase. Julia’s clearly integrated testing infrastructure has made our
development workflow significantly more robust.  Unicode support means that code
can correspond more closely to actual model equations, reducing the headache
associated with translating from "math" to "code".  Operator overloading and
user-defined syntax makes it easy to be much more expressive with our code. For
example, we can use `model[:α]` to access the value of parameter α from model
`model`, and `m <= Foo`  to add settings or parameters to the model object.
Inline markdown documentation also helps improve the developer experience.
Finally, we have found that Julia's highly integrated, Git-based package manager
provides a huge improvement over MATLAB's decentralized *FileExchange*. As Julia
users, we can now pull in high-quality, fully tested, community-supported
external packages that can each be installed or updated with a single command.
This reduces the incentive of users to create their own, lower-quality,
functionality, increasing developer *and* code performance.

We acknowledge that our package is far from perfect. Possible improvements to
*DSGE.jl* are many and varied. We may consider experimenting with alternative,
modern, numerical routines to improve speed. Ultimately, powerful
metaprogramming support (such as a macro to parse equilibrium conditions) would
allow user to specify model equations almost literally. We
[welcome](https://github.com/FRBNY-DSGE/DSGE.jl/blob/master/CONTRIBUTING.md)
such improvements to the existing code from the community.

## Challenges
Converting the FRBNY DSGE model from MATLAB, a mature and well-supported
language, to an extremely young language like Julia involved no shortage of
challenges. Significant changes to the Julia language itself are introduced in
rapid succession, and using *DSGE.jl* with a new Julia version inevitably floods
the user’s screen with deprecation warnings. There is significant difficulty in
finding written resources on the language in addition to the Julia Manual.
Google searches frequently return discussions in GitHub *Issues*, which are
unhelpful to elementary users and can be actively misleading at times.

Differences between the behavior of MATLAB and Julia’s core linear algebra
libraries led to many roadblocks in the development of *DSGE.jl*. Julia uses
multithreaded BLAS functions for some linear algebra functions.  Using a
different number of threads can change the results of matrix decomposition when
the matrix is singular. This indeterminacy caused significant problems for our
testing suite, both in comparing output matrices to MATLAB results and in
testing for reproducibility among Julia outputs.

We ran into many of these problems while porting the model solution algorithm,
`gensys`. At one point, the generalized Schur (QZ) decomposition is computed,
yielding the decompositions `A=QSZ'` and `B=QTZ'`. In MATLAB, upper triangular
matrices `S` and `T` are returned. In Julia, meanwhile, the default behavior is
to return a real decomposition with upper Hessenberg (blocked diagonal)
matrices `S` and `T`. [Add ordered QZ stuff.]

Finally, dealing with a recently introduced language can make it more
difficult for new users to produce performant code.  A typical economist,
especially one coming from a MATLAB background, may be unfamiliar with the
nature and use of language concepts like type stability, parametric types, and
preallocation. Julia’s profiler and debugger lack the flexibility of those in
MATLAB, and can make it difficult to identify the source of errors or
performance bottlenecks. However, Julia syntax should be familiar and intuitive
enough to be picked up quickly by MATLAB programmers without too much immediate
overhead.

## Conclusion
TODO

## Disclaimer
This post reflects the experience of the authors with Julia and MATLAB and does
not represent an endorsement by the Federal Reserve Bank of New York or the
Federal Reserve System of any particular product or service. The views
expressed in this post are those of the authors and do not necessarily reflect
the position of the Federal Reserve Bank of New York or the Federal Reserve
System. Any errors or omissions are the responsibility of the authors.
