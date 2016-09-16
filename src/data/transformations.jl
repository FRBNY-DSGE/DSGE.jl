"""
`annualtoquarter(v)`

Convert from annual to quarter frequency... by dividing by 4.
"""
function annualtoquarter(v)
    v / 4
end

"""
`quartertoannual(v)`

Convert from quarter to annual frequency... by multiplying by 4.
"""
function quartertoannual(v)
    4 * v
end


"""
`nominal_to_real(col, df; deflator_mnemonic=:GDPCTPI)`

Converts nominal to real values using the specified deflator.

## Arguments
- `col`: symbol indicating which column of `df` to transform
- `df`: DataFrame containining series for proper population measure and `col`

## Keyword arguments

- `deflator_mnemonic`: indicates which deflator to use to calculate real values. Default
  value is the FRED GDP Deflator mnemonic.
"""
function nominal_to_real(col, df; deflator_mnemonic=:GDPCTPI)
    return df[col] ./ df[deflator_mnemonic]
end


"""
```
percapita(m, col, df)
percapita(col, df, population_mnemonic)
```

Converts data column `col` of DataFrame `df` to a per-capita value.

## Arguments
- `col`: symbol indicating which column of data to transform
- `df`: DataFrame containining series for proper population measure and `col`
- `population_mnemonic`: a mnemonic found in df for some population measure.
"""
function percapita(m::AbstractModel, col::Symbol, df::DataFrame)
    population_mnemonic = get_setting(m, :population_mnemonic)
    percapita(col, df, population_mnemonic)
end
function percapita(col::Symbol, df::DataFrame, population_mnemonic::Symbol)
    df[col] ./ df[population_mnemonic]
end

"""
```
yt, yf = hpfilter(y, λ::Real)
```

Applies the Hodrick-Prescott filter ("H-P filter"). The smoothing parameter `λ` is applied
to the columns of `y`, returning the trend component `yt` and the cyclical component `yf`.  
For quarterly data, one can use λ=1600.

Consecutive missing values at the beginning or end of the time series are excluded from the
filtering. If there are missing values within the series, the filtered values are all NaN.

See also:
```
Hodrick, Robert; Prescott, Edward C. (1997). "Postwar U.S. Business Cycles: An Empirical
Investigation". Journal of Money, Credit, and Banking 29 (1): 1–16.
```
"""
function hpfilter(y, λ::Real)
    # Convert y to vector
    if !isa(y, Vector)
        try
            y = vec(y)
        catch
            error("Series must be convertible to Vector")
        end
    end

    # Indices of consecutive NaN elements at beginning
    i = 1
    j = length(y)
    while isnan(y[i])
        i = i+1
    end
    while isnan(y[j])
        j = j-1
    end

    # Filter and adjust for NaNs
    yt_, yf_ = hpfilter_(y[i:j], λ)
    yt = [fill(NaN, i-1); yt_; fill(NaN, length(y)-j)]
    yf = [fill(NaN, i-1); yf_; fill(NaN, length(y)-j)]

    return yt, yf
end

function hpfilter_{T<:Real}(y::Vector{T}, λ::Real)
    n = length(y)
    a = spzeros(n,n)
    for i = 3:n-2
        a[i,i]   = 6λ+1
        a[i,i+1] = -4λ
        a[i,i+2] = λ
        a[i,i-2] = λ
        a[i,i-1] = -4λ
    end

    a[2,2] = 1+5λ
    a[2,3] = -4λ
    a[2,4] = λ
    a[2,1] = -2λ
    a[1,1] = 1+λ
    a[1,2] = -2λ
    a[1,3] = λ

    a[n-1,n-1] = 1+5λ
    a[n-1,n-2] = -4λ
    a[n-1,n-3] = λ
    a[n-1,n]   = -2λ
    a[n,n]     = 1+λ
    a[n,n-1]   = -2λ
    a[n,n-2]   = λ

    yt = a\y
    yf = y-yt

    return yt, yf
end

"""
```
difflog(x::Vector{AbstractFloat})
```
"""
function difflog{T<:AbstractFloat}(x::Vector{T})
    [NaN; log(x[2:end]) - log(x[1:end-1])]
end


"""
```
difflog(x::DataArray{AbstractFloat})
```
"""
function difflog(x::DataArray)
    DSGE.na2nan!(x)
    y = convert(Vector{Float64}, x)
    return difflog(y)
end


"""
```
oneqtrpctchange(y)
```

Calculates the quarter-to-quarter percentage change of a series.
"""
function oneqtrpctchange(y)
    100 * difflog(y)
end


"""
```
hpadjust(y, df)
```

Adjust series to compensate for differences between filtered and unfiltered population.
## Arguments
- `y`: A vector of data
- `df`: DataFrame containing both a filtered and unfiltered population growth series
"""
function hpadjust(y, df; filtered_mnemonic=:filtered_population_growth,
                         unfiltered_mnemonic=:unfiltered_population_growth)
    y + 100 * (df[unfiltered_mnemonic] - df[filtered_mnemonic])
end
