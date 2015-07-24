# log posterior = log likelihood + log prior
# log Pr(Θ|YY)  = log Pr(YY|Θ)   + log Pr(Θ)
function posterior{T<:FloatingPoint}(model::AbstractModel, YY::Array{T, 2})
    return likelihood(model, YY) + prior(model.Θ)
end
    


# This is a dsge likelihood function that can handle 2-part estimation where
# there is a model switch.
# If there is no model switch, then we filter over the main sample all at once.
function likelihood{T<:FloatingPoint}(model::AbstractModel, YY::Array{T, 2})
    
    spec = model.spec



    # TODO: Find a way to return these matrices from measurement equation or something so `solve` isn't called twice
    ## step 1: solution to DSGE model - delivers transition equation for the state variables  S_t
    ## transition equation: S_t = TC+TTT S_{t-1} +RRR eps_t, where var(eps_t) = QQ
    TTT, CCC, RRR = solve(model::AbstractModel)


    
    ## step 2: define the measurement equation: X_t = ZZ*S_t + D + u_t
    ## where u_t = eta_t+MM* eps_t with var(eta_t) = EE
    ## where var(u_t) = HH = EE+MM QQ MM', cov(eps_t,u_t) = VV = QQ*MM'

    # Get measurement equation matrices set up
    #try
        ZZ, DD, QQ, EE, MM = measurement(model)
    #catch 
        # Error thrown during gensys
    #    return -Inf
    #end

    HH = EE + MM*QQ*MM'
    VV = QQ*MM'
    VVall =  [[RRR*QQ*RRR' RRR*VV]; [VV'*RRR' HH]]

    #=
    # TODO: What's happening here? In the Matlab code, this was only done for period 2
    if any(CCC != 0)
        DD = DD + ZZ*((UniformScaling(1) - TTT)\CCC)
    end
    =#


    
    ## step 3: compute log-likelihood using Kalman filter - written by Iskander
    ##         note that Iskander's program assumes a transition equation written as:
    ##         S_t = TTT S_{t-1} +eps2_t, where eps2_t = RRReps_t
    ##         therefore redefine QQ2 = var(eps2_t) = RRR*QQ*RRR'
    ##         and  VV2 = cov(eps2_t,u_u) = RRR*VV
    ##         define VVall as the joint variance of the two shocks VVall = var([eps2_tu_t])

    # TODO: Can we solve Lyapunov equation on matrices with anticipated shocks?
    # Solve lyapunov with normal period state matrices (i.e. period 2 matrices)
    A0 = zeros(spec["n_states_aug"], 1)
    P0 = dlyap!(copy(TTT), copy(RRR*QQ*RRR'))

    pyt, zend, Pend = kalcvf2NaN(YY', 1, zeros(spec["n_states_aug"], 1), TTT, DD, ZZ, VVall, A0, P0, 2)
    return pyt
end





# DLYAP   Discrete Lyapunov equation solver.
#
#    X = DLYAP(A,Q) solves the discrete Lyapunov equation:
#
#     A*X*A' - X + Q = 0
#
#    See also  LYAP.

#  J.N. Little 2-1-86, AFP 7-28-94
#  Copyright 1986-2001 The MathWorks, Inc.
#  $Revision: 1.11 $  $Date: 2001/01/18 19:50:01 $

#LYAP  Solve continuous-time Lyapunov equations.
#
#   x = LYAP(a,c) solves the special form of the Lyapunov matrix
#   equation:
#
#           a*x + x*a' = -c
#
#   See also  DLYAP.

#  S.N. Bangert 1-10-86
#  Copyright 1986-2001 The MathWorks, Inc.
#  $Revision: 1.10 $  $Date: 2001/01/18 19:50:23 $
#  Last revised JNL 3-24-88, AFP 9-3-95

# How to prove the following conversion is true.  Re: show that if
#         (1) Ad X Ad' + Cd = X             Discrete lyaponuv eqn
#         (2) Ac = inv(Ad + I) (Ad - I)     From dlyap
#         (3) Cc = (I - Ac) Cd (I - Ac')/2  From dlyap
# Then
#         (4) Ac X + X Ac' + Cc = 0         Continuous lyapunov
#
# Step 1) Substitute (2) into (3)
#         Use identity 2*inv(M+I) = I - inv(M+I)*(M-I)
#                                 = I - (M-I)*inv(M+I) to show
#         (5) Cc = 4*inv(Ad + I)*Cd*inv(Ad' + I)
# Step 2) Substitute (2) and (5) into (4)
# Step 3) Replace (Ad - I) with (Ad + I -2I)
#         Replace (Ad' - I) with (Ad' + I -2I)
# Step 4) Multiply through and simplify to get
#         X -inv(Ad+I)*X -X*inv(Ad'+I) +inv(Ad+I)*Cd*inv(Ad'+I) = 0
# Step 5) Left multiply by (Ad + I) and right multiply by (Ad' + I)
# Step 6) Simplify to (1)
function dlyap!(a, c)
    m, n = size(a)
    a = (a + UniformScaling(1))\(a - UniformScaling(1))
    c = (UniformScaling(1)-a)*c*(UniformScaling(1)-a')/2

    mc, nc = size(c)

    # a and c must be square and the same size
    if (m != n) || (m != mc) || (n != nc)
        error("Dimensions do not agree.")
    elseif m == 0
        x = zeros(m, m)
        return x
    end

    # Perform schur decomposition on a (and convert to complex form)
    ta, ua, _ = schur(complex(a)) # matlab code converts to complex - come back to
    # ua, ta = rsf2csf(ua, ta)
    # Schur decomposition of a' can be calculated from that of a.
    j = m:-1:1
    ub = ua[:, j]
    tb = ta[j ,j]'
    
    # Check all combinations of ta(i, i)+tb(j, j) for zero
    p1 = diag(ta).' # Use .' instead of ' in case a and a' are not real
    p2 = diag(tb)
    p_sum = abs(p1) .+ abs(p2)
    if any(p_sum .== 0) || any(abs(p1 .+ p2) .< 1000*eps()*p_sum)
        error("Solution does not exist or is not unique.")
    end

    # Transform c
    ucu = -ua'*c*ub

    # Solve for first column of transformed solution
    y = complex(zeros(n, n))
    y[:, 1] = (ta + UniformScaling(tb[1, 1]))\ucu[:, 1]

    # Solve for remaining columns of transformed solution
    for k=1:n-1
        km1 = 1:k
        y[:, k+1] = (ta + UniformScaling(tb[k+1, k+1]))\(ucu[:, k+1] - y[:, km1]*tb[km1, k+1])
    end

    # Find untransformed solution
    x = ua*y*ub'

    # Ignore complex part if real inputs (better be small)
    if isreal(a) && isreal(c)
        x = real(x)
    end

    # Force x to be symmetric if c is symmetric
    if issym(c)
        x = (x+x')/2
    end

    return x
end
