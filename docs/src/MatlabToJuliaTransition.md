# The DSGE MATLAB to Julia Transition: Improvements and Challenges

*Zac Cranko, Pearl Li, Spencer Lyon, Erica Moszkowski, Micah Smith, Pablo Winant*  

*December 3, 2015*

The FRBNY DSGE model is a relatively large New Keynesian model augmented with
financial frictions and a variety of innovations. Here at the Fed, we use it
both for forecasting and policy analysis. Research using this model includes
looking at the dynamics of [inflation during the great recession](https://www.newyorkfed.org/medialibrary/media/research/staff_reports/sr618.pdf),
the effects of [forward guidance](https://www.newyorkfed.org/medialibrary/media/research/staff_reports/sr574.pdf),
[and](http://libertystreeteconomics.newyorkfed.org/2015/05/why-are-interest-rates-so-low.html)
[much](http://libertystreeteconomics.newyorkfed.org/2014/09/an-assessment-of-the-frbny-dsge-models-real-time-forecasts.html)
[more](https://www.newyorkfed.org/medialibrary/media/research/staff_reports/sr695.pdf).

When we were approached by the folks at
[QuantEcon](http://quantecon.org/)
about a possible collaboration, we jumped at the idea, as it would give us an
opportunity to rework our code in an arguably faster language, redesign it from
the ground up, and release it open source for the benefit of the community. A
full-fledged package for the FRBNY DSGE model would also provide QuantEcon
another opportunity to highlight its contribution to high-performance,
quantitative economic modeling. Julia was the language of choice, recommended by
the QuantEcon group for its high performance and suitability for this breed of
technical computing.

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
to 3/4th that of the MATLAB code. The Metropolis-Hastings sampling step is the
most time consuming, and hence the relevant one in terms of assessing speed
improvement. On the basis of this step, we conclude that *DSGE.jl* is
approximately ten times faster than the MATLAB code.

How much of this increase is due to native performance adventures of Julia, and
how much is simply due to the improvements in design that came from rebuilding
this project from the ground up? It is of course difficult to say, and it is
important to emphasize that one cannot be sure what portion of the performance
increase can be attributed to inherent language features as opposed to design
differences. Indeed, our MATLAB code suffers from many inefficiencies due to its
long, cumulative development, and support for a plethora of models and features.
Meanwhile, these design issues have been largely addressed in our Julia package.
To best isolate differences in the languages themselves, we can look at our code
to compute the model solution with `gensys` and apply the Kalman filter with
`kalman_filter`. These two functions have relatively little redesign and
optimization as compared to the MATLAB code and provide the most comparable,
though still imperfect, measurements of performance. The reduction of 1/5th to
3/4th in computing time, therefore, could be taken as a first estimate of
Julia's advantage in this single arena of computation.

## Code Improvements

Julia provides versatile language features that allow us to improve our code's
performance and clarity in several fundamental ways. First and foremost of these
is the highly integrated, robust, and flexible type system that lends itself
naturally to our DSGE model. At the center of the *DSGE.jl* package is the
*model object*. Here, one can store all information associated with the model –
including the numerous parameters, priors, states, equilibrium conditions,
computational settings, and flags – in one place.  By simply passing the model
object as an argument to any function, the function has access to all of the
model's fields.  By comparison, our MATLAB code stored all variables
directly in the global workspace – an approach that scaled poorly as model
specifications become more and more complex. To illustrate just how unwieldy our
MATLAB code was, many of our function calls required more than 20 positional
arguments, a serious challenge for usage and human-readability:
```matlab
[post_new,like_new,zend_new,ZZ_new,DD_new,QQ_new] = ...
    feval('objfcnmhdsge',para_new,bounds,YY,YY0,nobs,nlags,nvar,mspec,npara,...
    trspec,pmean,pstdd,pshape,TTT_new,RRR_new,CCC_new,valid_new,para_mask,...
    coint,cointadd,cointall,YYcoint0,args_nant_antlags{:});
```
While several of these arguments (e.g., `coint`) relate to a feature not-implemented
in Julia, one can still see the excesses of providing so much information about
the model separately in function calls.

The same code in Julia:
```julia
post_out = posterior!(m, para_new, data; mh=true)
```

Certainly, one could approximate a "model object" in MATLAB by using its own
object-oriented classes, or by "bundling" model attributes into a `struct` or
other data structure.  However, MATLAB classes are both relatively complicated
and slower than non-object implementations. And using `struct`s in this way
results in copies of all model variables made on every function call.

Indeed, changes like this reduce the number of lines of code in *DSGE.jl*, a
rough proxy for ease of maintenance. We find that the fixed cost of setting up
the type system is offset by savings in core programs.

| Language             | Lines of code (hundreds) |
| -------------------- | :----------------------: |
| Matlab               | 63                       |
| Julia                | 37                       |

A type-based approach allows us to take advantage of method dispatch in Julia by
defining different model types for different model specifications. As detailed
in the
[README file](https://github.com/FRBNY-DSGE/DSGE.jl/blob/master/README.md),
changes to the model's equilibrium conditions and measurement equation are
referred to as changes in a model's "specification."  In the Julia code, model
specifications have a 1:1 correspondence with concrete types.  Where necessary,
a single *function* can have multiple *methods* defined, that are customized for
different model types. For example, the `augment_states` function
augments the model's transition matrices after it has been solved.  We can
pass any model object `m` to `augment_states`, and Julia ensures that the proper,
model-specific *method* is *dispatched*:
```julia
TTT, RRR, CCC = augment_states(m, TTT_gensys, RRR_gensys, CCC_gensys)
```

In our MATLAB code, on the other hand, we would approximate this type of
dispatch by using a `switch` statement over a model identifier variable. For
the hundreds of models we have worked with in a development capacity, this led
to bloat in our model solution code. In Julia, we encapsulate this behavior
within the model definition itself.

It is easy to see that all model types constructed for use with *DSGE.jl* are
closely related: they will have the same fields, and are passed to the same
methods.  If it sounds to you like we have an implicit interface here, you’re
right. Rather than implementing each object as a standalone type, we define an
abstract type, `AbstractModel`, to serve as the parent for all model types.
Because most of our routines are not model-specific, we need only define them
once (with an argument of type `AbstractModel`) and Julia's dispatch system
takes care of the rest. We similarly define model parameters as subtypes of a
common abstract type `AbstractParameter`. This allows us to abstract to one
notion of a model parameter, while implementing different kinds of parameters in
different ways. We also use parameterized (generic) types to increase the
flexibility of model parameters (as well as elsewhere in our codebase):
```julia
# A parameter contains values of data type `T` and embeds a transformation of
# type `U`
abstract Parameter{T,U<:Transform} <: AbstractParameter{T}

# One transformation used is the identity
type UnscaledParameter{T,U} <: Parameter{T,U}
    # ...
end
```

These functions expect the model object to have certain fields,
and for those fields to have certain types. (As an example of Julia's
youthful status as a language,
[discussion continues](https://github.com/JuliaLang/julia/issues/6975),
as of this writing, on an appropriate manner to explicitly introduce
interfaces.)

With a clear interface in place, running new model specifications
using *DSGE.jl* is relatively straightforward. (See
[here](@ref editing-extending-model)
for detailed instructions).

Julia's JIT compilation provides significant performance boosts in some
areas. For example, we allow a variable number of anticipated monetary policy
shocks, beginning in 2008Q4, that we use to treat the zero lower bound. In our
MATLAB code, we suffer some dynamic code generation to implement this feature.
```matlab
if exist('nant','var')
    for i=1:nant
        eval(strcat('rm_tl',num2str(i),'  = ',num2str(nstates+i)));
        eval(strcat('rm_tl',num2str(i),'  = ',num2str(nstates+i)));
    end
end
```

Julia's faster evaluation of such statements reduces this performance
hit, as these symbols can be associated with the model object.
```julia
[symbol("rm_tl$i") for i = 1:n_anticipated_shocks(m)]
# ...
[symbol("rm_shl$i") for i = 1:n_anticipated_shocks(m)]
```

Granted, there may be better solutions to our problem in both languages,
but similar situations involving code generation are easily addressed in Julia.

We have found that a number of Julia features make working with
*DSGE.jl* simply more pleasant and user-friendly than working with our old
codebase. Julia's clearly integrated testing infrastructure has made our
development workflow significantly more robust. Unicode support means that code
can correspond more closely to actual model equations, reducing the headache
associated with translating from "math" to "code".  (Inline
[Markdown](https://en.wikipedia.org/wiki/Markdown)
documentation helps in a similar way.) Operator overloading and user-defined
syntax make it easy to be much more expressive with our code.
```julia
julia> m[:α]                         # Access value of param α from model m
julia> m <= parameter(:ϵ_p, 10.000)  # Add parameter ϵ_p to model
julia> Γ0, Γ1, C, Ψ, Π  = eqcond(m)  # Get equilibrium conditions
```

We have found that Julia's highly integrated, Git-based package manager
is an improvement over MATLAB's decentralized
[FileExchange](http://www.mathworks.com/matlabcentral/fileexchange/).
As Julia users, we can now pull in high-quality, fully tested,
community-supported external packages that can each be installed or updated with
a single command.
```julia
julia> Pkg.add("QuantEcon")          # That's it!
```
This reduces the need for users to create their own, likely lower-quality
functionality, increasing developer *and* code performance. (Or the need to
fight for the toolbox licenses available to their department.)

We acknowledge that our package is far from perfect. Possible improvements to
*DSGE.jl* are many and varied. We may experiment with alternative, modern,
numerical routines to improve speed. Ultimately, powerful metaprogramming
support would allow user to specify model equations more literally, in
mathematical notation. We
[welcome](@ref contributing) improvements to the existing code from the community.

## Challenges

Converting the FRBNY DSGE model from MATLAB, a mature and well-supported
language, to an extremely young language like Julia involved no shortage of
challenges. Significant changes to the Julia language itself are introduced in
rapid succession, and using *DSGE.jl* with a new Julia version inevitably floods
the user’s screen with deprecation warnings. There is significant difficulty in
finding written resources on the language beyond the Julia Manual itself.
Google searches frequently return discussions in GitHub *Issues*, which are
unhelpful to elementary users and can be actively misleading at times.

Differences between the behavior of MATLAB and Julia’s core linear algebra
libraries led to many roadblocks in the development of *DSGE.jl*. Julia uses
multithreaded BLAS functions for some linear algebra functions.  Using a
different number of threads can change the results of matrix decomposition when
the matrix is singular. This indeterminacy caused significant problems for our
testing suite, both in comparing output matrices to MATLAB results and in
testing for reproducibility among Julia outputs.

We ran into similar numerical problems while porting the model solution
algorithm, `gensys`. At one point, the generalized Schur (QZ) decomposition is
computed, yielding the decompositions `A=QSZ'` and `B=QTZ'`. In MATLAB, upper
triangular matrices `S` and `T` are returned. In Julia, meanwhile, the default
behavior is to return a real decomposition with upper Hessenberg (blocked
diagonal) matrices `S` and `T`. Differing behaviors like this in the two
languages might expose a user without deep knowledge of the procedure to
errors.

Finally, dealing with a recently introduced language can make it more
difficult for new users to produce performant code.  A typical economist,
especially one coming from a MATLAB background, may be unfamiliar with the
nature and use of language concepts like type stability, parametric types, and
preallocation. Julia's profiler and debugger lack the flexibility of those in
MATLAB, and can make it difficult to identify the source of errors or
performance bottlenecks. And Julia IDEs, like [Juno](http://junolab.org/), while
admirable, are not as mature or featured as the MATLAB IDE.

It is important to note again that similar improvements could have been made to
our MATLAB code directly. (As we would be the first ones to admit.) Regardless,
the Julia paradigm results in code that is high-quality from the outset.

## Conclusion

After months of hard work, we are pleased to be able to increase the performance
of our model and provide our project for the benefit of the community. For those
considering similar projects, we find the benefits of a transition to Julia are
significant. One should, however, be realistic about the challenges that will be
faced transitioning to a young language.

## Disclaimer

This post reflects the experience of the authors with Julia and MATLAB and does
not represent an endorsement by the Federal Reserve Bank of New York or the
Federal Reserve System of any particular product or service. The views
expressed in this post are those of the authors and do not necessarily reflect
the position of the Federal Reserve Bank of New York or the Federal Reserve
System. Any errors or omissions are the responsibility of the authors.
