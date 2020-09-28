#=
This code is based on a routine originally copyright Chris Sims.
See http://sims.princeton.edu/yftp/optimize/
=#

struct Csminwel <: Optim.SecondOrderOptimizer end

const rc_messages = Dict(0 => "Standard Iteration",
                         1 => "zero gradient",
                         2 => "back and forth on step length never finished",
                         3 => "smallest step still improving too slow",
                         4 => "back and forth on step length never finished",
                         5 => "largest step still improving too fast",
                         6 => "smallest step still improving too slow, reversed gradient",
                         7 => "warning: possible inaccuracy in H matrix")

macro csminwelltrace()
    esc(quote
        if tracing
            dt = Dict()
            if extended_trace
                dt["x"] = copy(x)
                dt["g(x)"] = copy(gr)
                dt["H(x)"] = copy(H)
                dt["rc"] = copy(rc)
                dt["rc_message"] = rc_messages[rc]
            end
            grnorm = norm(gr, Inf)
            Optim.update!(tr,
                          iteration,
                          f_x,
                          grnorm,
                          dt,
                          store_trace,
                          show_trace)
        end
    end)
end

"""
```
csminwel(fcn::Function, grad::Function, x0::Vector, H0::Matrix=1e-5.*eye(length(x0)), args...;
         xtol::Real=1e-32, ftol::Float64=1e-14, grtol::Real=1e-8, iterations::Int=1000,
         store_trace::Bool = false, show_trace::Bool = false, extended_trace::Bool = false,
         verbose::Symbol = :none, rng::AbstractRNG = MersenneTwister(0), kwargs...)
```

Minimizes `fcn` using the csminwel algorithm.

### Arguments

* `fcn::Function` : The objective function
* `grad::Function` : The gradient of the objective function. This argument can be omitted if
an analytical gradient is not available, which will cause a numerical gradient to be
calculated.
* `x0::Vector`: The starting guess for the optimizer

### Optional Arguments

* `H0::Matrix`: An initial guess for the Hessian matrix -- must be
positive definite. If none is given, then a scaled down identity
matrix is used.
* `args...`:  Other positional arguments to be passed to `f` on each
function call

### Keyword Arguments

* `ftol::{T<:Real}=1e-14`: Threshold for convergence in terms of change
in function value across iterations.
* `iterations::Int=100`: Maximum number of iterations
* `kwargs...`: Other keyword arguments to be passed to `f` on each
function call
"""
function csminwel(fcn::Function,
                  grad::Function,
                  x0::Vector,
                  H0::Matrix=1e-5.*eye(length(x0)),
                  args...;
                  xtol::Real           = 1e-32,  # default from Optim.jl
                  ftol::Float64        = 1e-14,  # Default from csminwel
                  grtol::Real          = 1e-8,   # default from Optim.jl
                  iterations::Int      = 1000,
                  store_trace::Bool    = false,
                  show_trace::Bool     = false,
                  extended_trace::Bool = false,
                  verbose::Symbol      = :none,
                  rng::AbstractRNG     = MersenneTwister(0),
                  kwargs...)

    if show_trace
        @printf "Iter     Function value   Gradient norm \n"
    end

    # unpack dimensions
    nx = size(x0, 1)

    # Count function and gradient calls
    f_calls = 0
    g_calls = 0

    # Maintain current state in x and previous state in x_previous
    x, x_previous = copy(x0), copy(x0)

    # start with Initial Hessian
    H = H0

    # start rc parameter at 0
    rc = 0

    f_x = fcn(x0, args...; kwargs...)
    f_calls += 1

    if isinf(f_x)
        throw(ArgumentError("Bad initial guess, `f_x` returned Inf. Try again"))
    elseif f_x > 1e50
        throw(ArgumentError("Bad initial guess. Try again"))
    end

    gr, badg = grad(x0)
    g_calls += 1

    # Count iterations
    iteration = 0

    # Maintain a trace
    tr = OptimizationTrace{Float64, Csminwel}()
    tracing = show_trace || store_trace || extended_trace
    @csminwelltrace

    # set objects to their starting values
    retcode3 = 101

    # set up return variables so they are available outside while loop
    fh = copy(f_x)
    xh = copy(x0)
    gh = copy(x0)
    retcodeh = 1000

    # Assess multiple types of convergence
    x_converged, f_converged, gr_converged = false, false, false

    # Declare residual variables
    x_resid, f_resid, gr_resid = 0.0, 0.0, 0.0

    # Iterate until convergence or exhaustion
    converged = false
    while !converged && iteration < iterations
        iteration += 1

        f1, x1, fc, retcode1 = csminit(fcn, x, f_x, gr, badg, H, args...;
                                       verbose=verbose, kwargs...)
        f_calls += fc

        if retcode1 != 1
            if retcode1 == 2 || retcode1 == 4
                wall1, badg1 = true, true
            else
                g1, badg1 = grad(x1)
                g_calls += 1
                wall1 = badg1
            end

            # Bad gradient or back and forth on step length.  Possibly at cliff edge. Try
            # perturbing search direction if problem not 1D
            if wall1 && (length(H) > 1)
                println(verbose, :low, "Cliff. Perturbing search direction.")

                Hcliff = H + Matrix(Diagonal(diag(H).*rand(rng, nx)))
                f2, x2, fc, retcode2 = csminit(fcn, x, f_x, gr, badg, Hcliff,
                                               args...; verbose=verbose, kwargs...)
                f_calls += fc

                if f2 < f_x
                    if retcode2==2 || retcode2==4
                        wall2 = true; badg2 = true
                    else
                        g2, badg2 = grad(x2)
                        g_calls += 1
                        wall2 = badg2
                        badg2
                    end

                    if wall2
                        println(verbose, :low, "Cliff again. Try traversing")

                        if norm(x2-x1) < 1e-13
                            f3 = f_x
                            x3 = x
                            badg3 = true
                            retcode3 = 101
                        else
                            gcliff = ( (f2-f1) / ((norm(x2-x1))^2) )*(x2-x1)
                            if (size(x0 , 2)>1)
                                gcliff = gcliff'
                            end
                            f3, x3, fc, retcode3 = csminit(fcn, x, f_x, gcliff,
                                                           false, eye(nx),
                                                           args...; verbose=verbose, kwargs...)
                            f_calls += fc

                            if retcode3==2 || retcode3==4
                                wall3 = true
                                badg3 = true
                            else
                                g3, badg3 = grad(x3)
                                g_calls += 1
                                wall3 = badg3
                            end
                        end
                    else
                        f3 = f_x
                        x3 = x
                        badg3 = true
                        retcode3 = 101
                    end
                else
                    f3 = f_x
                    x3 = x
                    badg3 = true
                    retcode3 = 101
                end
            else
                # normal iteration, no walls, or else 1D, or else we're finished here.
                f2, f3 = f_x, f_x
                badg2, badg3 = true, true
                retcode2, retcode3 = 101, 101
            end
        else
            f1, f2, f3 = f_x, f_x, f_x
            retcode2, retcode3 = retcode1, retcode1
        end

        # how to pick gh and xh
        if f3 < f_x - ftol && badg3==0
            ih = 3
            fh = f3
            xh = x3
            gh = g3
            badgh = badg3
            retcodeh = retcode3
        elseif f2 < f_x - ftol && badg2==0
            ih = 2
            fh = f2
            xh = x2
            gh = g2
            badgh = badg2
            retcodeh = retcode2
        elseif f1 < f_x - ftol && badg1==0
            ih = 1
            fh = f1
            xh = x1
            gh = g1
            badgh = badg1
            retcodeh = retcode1
        else
            fh, ih = findmin([f1 , f2 , f3])

            if ih == 1
                xh = x1
                retcodeh = retcode1
            elseif ih == 2
                xh = x2
                retcodeh = retcode2
            elseif ih == 3
                xh = x3
                retcodeh = retcode3
            end

            if @isdefined gh
                nogh = isempty(gh)
            else
                nogh = true
            end

            if nogh
                gh, badgh = grad(xh)
                g_calls += 1
            end

            badgh = true
        end

        stuck = (abs(fh-f_x) < ftol)
        if !badg && !badgh && !stuck
            H = bfgsi(H , gh-gr , xh-x; verbose=verbose)
        end

        println(verbose, :high, @sprintf "Improvement on iteration %d = %18.9f\n" iteration fh-f_x)

        if stuck
            println(verbose, :low, "improvement < ftol -- terminating")
        end

        # record# retcodeh of previous x
        copyto!(x_previous, x)

        # update before next iteration
        f_x_previous, f_x = f_x, fh
        x = xh
        gr = gh
        badg = badgh

        # Check convergence
        x_converged,
        f_converged,
        gr_converged,
        converged = assess_convergence(x,
                                       x_previous,
                                       f_x,
                                       f_x_previous,
                                       gr,
                                       xtol,
                                       ftol,
                                       grtol)

        x_resid  = Optim.x_abschange(x, x_previous)
        f_resid  = Optim.f_abschange(f_x, f_x_previous)/abs(f_x + ftol)
        gr_resid = Optim.g_residual(gr)

        @csminwelltrace
    end
    try
        # Optim v1.x adds an additional field to the MultivariateOptimizationResults type
        return MultivariateOptimizationResults(Csminwel(), x0, x, convert(Float64, f_x),
                                               iteration, iteration==iterations, x_converged, xtol, xtol, x_resid, x_resid,
                                               f_converged, ftol, ftol, f_resid, f_resid,
                                               gr_converged, grtol, gr_resid, false, tr, f_calls, g_calls, 0,
                                               false, NaN, NaN, nothing), H  # also return H
    catch e
        if isa(e, MethodError) # Then likely using an older version of Optim
            return MultivariateOptimizationResults(Csminwel(), x0, x, convert(Float64, f_x),
                                                   iteration, iteration==iterations, x_converged, xtol, xtol, x_resid, x_resid,
                                                   f_converged, ftol, ftol, f_resid, f_resid,
                                                   gr_converged, grtol, gr_resid, false, tr, f_calls, g_calls, 0,
                                                   false, NaN, NaN), H  # also return H
        else
            rethrow(e)
        end
    end

# the NaNs are for time and timelimit--apparently if they're NaN, there is no time limit
# the last false if ls_success=false...not quite sure?
    #=return MultivariateOptimizationResults(Csminwel(), x0, x, convert(Float64, f_x),
        iteration, iteration==iterations, x_converged, xtol, x_resid, f_converged, ftol, f_resid,
        gr_converged, grtol, gr_resid, false, tr, f_calls, g_calls, 0), H  # also return H=#
end


#=
Version of `csminwel` that will use finite differencing methods to
approximate the gradient numerically. This is convenient for cases where
you cannot supply an analytical derivative, but it is not as robust as
using the true derivative.
=#
function csminwel(fcn::Function,
                  x0::Vector,
                  H0::Matrix=0.5.*eye(length(x0)),
                  args...;
                  xtol::Real           = 1e-32, # default from Optim.jl
                  ftol::Float64        = 1e-14, # Default from csminwel
                  grtol::Real          = 1e-8,  # default from Optim.jl
                  iterations::Int      = 1000,
                  store_trace::Bool    = false,
                  show_trace::Bool     = false,
                  extended_trace::Bool = false,
                  verbose::Symbol      = :none,
                  rng::AbstractRNG     = MersenneTwister(0),
                  kwargs...)

    grad(x::Array{T}) where {T <: Number} = csminwell_grad(fcn, x, args...; kwargs...)
    csminwel(fcn, grad, x0, H0, args...;
             xtol=xtol, ftol=ftol, grtol=grtol, iterations=iterations,
             store_trace=store_trace, show_trace=show_trace,
             extended_trace=extended_trace, verbose=verbose, rng=rng, kwargs...)
end


function csminwell_grad(fcn, x, args...; kwargs...)
    f(a) = fcn(a, args...; kwargs...)
    gr = Calculus.gradient(f, x)
    bad_grads = abs.(gr) .>= 1e15
    gr[bad_grads] .= 0.0
    return gr, any(bad_grads)
end

function csminit(fcn, x0, f0, g0, badg, H0, args...; verbose::Symbol=:none, kwargs...)
    angle = .005

    #(0<THETA<.5) THETA near .5 makes long line searches, possibly fewer iterations.
    theta = .3
    fchange = 1000
    minlamb = 1e-9
    mindfac = .01
    f_calls = 0
    lambda = 1.0
    xhat = x0
    f = f0
    fhat = f0
    gr = g0
    gnorm = norm(gr)

    if gnorm < 1e-12 && !badg
        # gradient convergence
        retcode = 1
        dxnorm = 0.0
    else
        # with badg true, we don't try to match rate of improvement to directional
        # derivative.  We're satisfied just to get some improvement in f.
        dx = vec(-H0*gr)
        dxnorm = norm(dx)

        if dxnorm > 1e12
            println(verbose, :low, "Near singular H problem.")
            dx = dx * fchange / dxnorm
        end

        dfhat = dot(dx, g0)

        if !badg
            # test for alignment of dx with gradient and fix if necessary
            a = -dfhat / (gnorm*dxnorm)

            if a < angle
                dx -= (angle*dxnorm/gnorm + dfhat/(gnorm*gnorm)) * gr
                dx *= dxnorm/norm(dx)
                dfhat = dot(dx, gr)

                println(verbose, :low, @sprintf "Correct for low angle %f\n" a)
            end
        end

        println(verbose, :high, @sprintf "Predicted Improvement: %18.9f\n" -dfhat/2)

        # Have OK dx, now adjust length of step (lambda) until min and max improvement rate
        # criteria are met.
        done = false
        fact = 3.0
        shrink = true
        lambda_min = 0.0
        lambda_max = Inf
        lambda_peak = 0.0
        f_peak = f0
        lambda_hat = 0.0

        while !done
            if size(x0, 2) > 1
                dxtest = x0 + dx' * lambda
            else
                dxtest = x0 + dx * lambda
            end

            f = fcn(dxtest, args...; kwargs...)

            println(verbose, :high, @sprintf "lambda = %10.5f; f = %20.7f\n" lambda f)

            if f < fhat
                fhat = f
                xhat = dxtest
                lambdahat = lambda
            end

            f_calls += 1
            shrink_signal = (!badg &&
                             (f0-f < maximum([-theta*dfhat*lambda 0]))) ||
                            (badg && (f0-f) < 0)

            grow_signal = !badg && (lambda > 0)  &&
                          (f0-f > -(1-theta)*dfhat*lambda)

            if shrink_signal && ((lambda > lambda_peak) || lambda < 0 )
                if (lambda > 0) && ((!shrink) || (lambda/fact <= lambda_peak))
                    shrink = true
                    fact = fact^.6
                    while lambda / fact <= lambda_peak
                        fact = fact^.6
                    end

                    if abs(fact - 1.0) < mindfac
                        if abs(lambda) < 4
                            retcode = 2
                        else
                            retcode = 7
                        end

                        done = true
                    end
                end

                if lambda < lambda_max && lambda > lambda_peak
                    lambda_max = lambda
                end

                lambda /= fact
                if abs(lambda) < minlamb
                    if (lambda > 0) && (f0 <= fhat)
                        lambda = -lambda*fact^6
                    else
                        if lambda < 0
                            retcode = 6
                        else
                            retcode = 3
                        end
                        done = true
                    end
                end


            elseif (grow_signal && lambda > 0) || (shrink_signal &&
                                    ((lambda <= lambda_peak) && (lambda > 0)))
                if shrink
                    shrink = false
                    fact = fact^.6
                    if abs(fact - 1) < mindfac
                        if abs(lambda) < 4
                            retcode = 4
                        else
                            retcode = 7
                        end
                        done = true
                    end
                end

                if f < f_peak && lambda > 0
                    f_peak = f
                    lambda_peak = lambda
                    if lambda_max <= lambda_peak
                        lambda_max = lambda_peak * fact^2
                    end
                end

                lambda *= fact
                if abs(lambda) > 1e20
                    retcode = 5
                    done = true
                end
            else
                done = true
                if fact < 1.2
                    retcode = 7
                else
                    retcode = 0
                end
            end
        end
    end

    println(verbose, :high, @sprintf "Norm of dx %10.5f\n" dxnorm)

    return fhat, xhat, f_calls, retcode
end

"""
```
bfgsi(H0, dg, dx)
```

### Arguments
- `H0`: hessian matrix
- `dg`: previous change in gradient
- `dx`: previous change in x
"""
function bfgsi(H0, dg, dx; verbose::Symbol = :none)
    if size(dg, 2) > 1
        dg = dg'
    end

    if size(dx, 2) > 1
        dx = dx'
    end

    Hdg = H0*dg
    dgdx = dot(dx, dg)

    H = H0
    if abs(dgdx) > 1e-12
        H += (dgdx.+(dg'*Hdg)).*(dx*dx')/(dgdx^2) - (Hdg*dx'.+dx*Hdg')/dgdx
    elseif norm(dg) < 1e-7
        # gradient is super small so don't worry updating right now
        # do nothing
    else
        @warn "bfgs update failed"

        println(verbose, :high, @sprintf "|dg| = %f, |dx| = %f\n" norm(dg) norm(dx))
        println(verbose, :high, @sprintf "dg'dx = %f\n" dgdx)
        println(verbose, :high, @sprintf "|H*dg| = %f\n" norm(Hdg))
    end

    return H
end


function assess_convergence(x::Array,
                            x_previous::Array,
                            f_x::Real,
                            f_x_previous::Real,
                            gr::Array,
                            xtol::Real,
                            ftol::Real,
                            grtol::Real)
    x_converged, f_converged, gr_converged = false, false, false

    if Optim.maxdiff(x, x_previous) < xtol
        x_converged = true
    end

    # Relative Tolerance
    # if abs(f_x - f_x_previous) / (abs(f_x) + ftol) < ftol || nextfloat(f_x) >= f_x_previous
    # Absolute Tolerance
    if abs(f_x - f_x_previous) < ftol
        f_converged = true
    end

    if norm(vec(gr), Inf) < grtol
        gr_converged = true
    end

    converged = x_converged || f_converged || gr_converged

    return x_converged, f_converged, gr_converged, converged
end

# function Base.show(io::IO, t::OptimizationState)
#     @printf io "%6d   %14.10e   %14.10e\n" t.iteration t.value t.gradnorm
#     if !isempty(t.metadata)
#         for (key, value) in t.metadata
#             @printf io " * %s: %s\n" key value
#         end
#     end
#     return
# end
