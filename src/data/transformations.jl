function Base.call(x::Function, df::DataFrame)
    x(df)
end


"""
`function thousandstomillions(v)`

Convert vector from thousands scale to millions scale
"""
function thousandstomillions(v)
    v / 1000
end


"""
`function annualtoquarter(v)`

Convert from annual to quarter by dividing by 4
"""
function annualtoquarter(v)
    v / 4
end

function nominal_to_real(col, df, deflator=:GDP)
    if deflator == :GDP
        return df[col] ./ df[:GDPCTPI] 
    elseif deflator == :PCE
        return df[col] ./ df[:PCEPILFE]
    else
        error("nominal_to_real isn't defined for this deflator")
    end
end

"""
`function nominal_to_realpercapita(col::Symbol, df)`

Convert column col from a nominal to a real, per-capita value
"""
function nominal_to_realpercapita(col, df, deflator=:GDP)
    nominal_to_real(col, df, deflator) ./ thousandstomillions(df[:CNP16OV])
end

"""
"""
function to_perpersonemployed(col, df; src=:nonfarm)
    if src == :nonfarm
        df[col] ./ df[:PRS85006013]
    elseif src == :civilian
        df[col] ./ df[:CE16OV]
    else
        error("to_perpersonemployed isn't defined for this deflator")
    end
end


"""
`function nominal_to_realperemployed(col::Symbol, df)`

Convert column col from a nominal to a real, per-person-employed value
"""
function nominal_to_realperemployed(col, df; src=:nonfarm)
    if src == :nonfarm
        to_perpersonemployed(col, df, src=src) ./ df[:PRS85006013]
    elseif src == :civilian
        to_perpersonemployed(col, df, src=src) ./ df[:CE16OV]
    end
end


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


