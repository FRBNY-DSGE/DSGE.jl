"""
`dlyap(a, c)`

Discrete Lyapunov equation solver. The discrete Lyapunov equation is given by 
```
A*X*A' - X + C = 0
```

Attribution
-----------
Adapted from `DLYAP`, J.N. Little, The MathWorks Inc., 2001-01-18
"""
function dlyap(a, c)
    return dlyap!(copy(a), copy(c))
end
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
