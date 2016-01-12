"""
```
yt, yf = hpfilter(y::Vector{AbstractFloat}, λ::AbstractFloat)
```

Applies the Hodrick-Prescott filter. The smoothing parameter `λ` is applied to the columns of
`y`, returning the trend component `yt` and the cyclical component `yf`.
"""
function hpfilter(y::Vector{AbstractFloat}, λ::AbstractFloat)
    n = length(y);
    a = spzeros(n,n);
    for i = 3:n-2
        a(i,i)   = 6λ+1;
        a(i,i+1) = -4λ;
        a(i,i+2) = λ;
        a(i,i-2) = λ;
        a(i,i-1) = -4λ;
    end

    a(2,2) = 1+5λ;
    a(2,3) = -4λ;
    a(2,4) = λ;
    a(2,1) = -2λ;
    a(1,1) = 1+λ;
    a(1,2) = -2λ;
    a(1,3) = λ;

    a(n-1,n-1) = 1+5λ;
    a(n-1,n-2) = -4λ;
    a(n-1,n-3) = λ;
    a(n-1,n)   = -2λ;
    a(n,n)     = 1+λ;
    a(n,n-1)   = -2λ;
    a(n,n-2)   = λ;

    yt = a\y;
    yf = y-yt;
end

"""
```
difflog(y::Vector{AbstractFloat})
```
"""
function difflog(x::Vector{AbstractFloat})
    return log([NaN; diff(x)])
end
