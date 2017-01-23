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
`quartertoannualpercent(v)`

Convert from quarter to annual frequency in percent... by multiplying by 400.
"""
function quartertoannualpercent(v)
    400 * v
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
    population_mnemonic = parse_population_mnemonic(m)[1]
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




## REVERSE TRANSFORMS

"""
```
logtopct_annualized_percapita(y, pop_growth, q_adj = 100)
```
Transform from log growth rates to % growth rates (annualized).

### Note

This should only be used for output, consumption, investment
and GDP deflator (inflation).

### Inputs

- `y`: The data we wish to transform to 4 quarter annualized percent change from
  1-quarter log-levels. `y` is either a vector of length `nperiods` or an
  `ndraws x `nperiods` matrix.

- `pop_growth::Vector`: The length `nperiods` vector of population growth rates.

- `q_adj`: Optional argument defaulting to 100.
"""
function logtopct_annualized_percapita(y::Array, pop_growth::Vector, q_adj = 100.)
    # `y` is either a vector of length `nperiods` or an
    # `ndraws` x `nperiods` matrix
    if ndims(y) == 1
        nperiods = length(y)
    else
        nperiods = size(y, 2)

        # Transpose `pop_growth` to a 1 x `nperiods` row vector so it can be
        # broadcasted to match the dimensions of `y`
        pop_growth = pop_growth'
    end

    @assert length(pop_growth) == nperiods

    100. * (exp(y/q_adj .+ pop_growth).^4 - 1.)
end

"""
```
logtopct_annualized(y, q_adj = 100)
```

Transform from log growth rates to total (not per-capita) % growth
rates (annualized).
"""
function logtopct_annualized(y, q_adj = 100.)
    100. * (exp(y/q_adj).^4 - 1.)
end

"""
```
loglevelto4qpct_annualized(y, y0)
```

Transform from log level to 4-quarter annualized percent change

### Note

This is usually applied to labor supply (hours worked per hour), and
probably shouldn't be used for any other observables.

### Inputs

- `y`: The data we wish to transform to 4 quarter annualized percent change from
  1-quarter log-levels. `y` is either a vector of length `nperiods` or an
  `ndraws x `nperiods` matrix.

- `y0`: The last data point in the history (of state or observable)
  corresponding to the `y` variable.  This is required to compute a percentage
  change for the first period.
"""
function loglevelto4qpct_annualized{T<:AbstractFloat}(y::Array, y0::T)
    # `y_t1` is an array of the same size as `y`, representing the previous
    # period observations for each draw
    if ndims(y) == 1
        y_t1 = vcat([y0], y)
    else
        ndraws = size(y, 1)
        y0s  = fill(y0, ndraws, 1)
        y_t1 = hcat(y0s, y[:, 1:end-1])
    end

    # Subtract log levels to get log growth rates, then take the exponential to
    # get growth rates
    100. * (exp(y./100. - y_t1./100.).^4 .- 1.)
end

"""
```
loglevelto4qpct_annualized_percapita(y, y0, pop_growth)
```

Transform from log level to 4-quarter annualized percent change, adjusting for
population growth.

### Note

This is usually applied to labor supply (hours worked per hour), and
probably shouldn't be used for any other observables.

### Inputs

- `y`: The data we wish to transform to 4 quarter annualized percent change from
  1-quarter log-levels. `y` is either a vector of length `nperiods` or an
  `ndraws x `nperiods` matrix.

- `y0`: The last data point in the history (of state or observable)
  corresponding to the `y` variable.  This is required to compute a percentage
  change for the first period.

- `pop_growth::Vector`: The length `nperiods` vector of population growth rates.
"""
function loglevelto4qpct_annualized_percapita{T<:AbstractFloat}(y::Array, y0::T, pop_growth::Vector)
    # `y_t1` is an array of the same size as `y`, representing the previous
    # period observations for each draw
    if ndims(y) == 1
        nperiods = length(y)
        y_t1 = vcat([y0], y)
    else
        (ndraws, nperiods) = size(y)
        y0s  = fill(y0, ndraws, 1)
        y_t1 = hcat(y0s, y[:, 1:end-1])

        # Transpose `pop_growth` to a 1 x `nperiods` row vector so it can be
        # broadcasted to match the dimensions of `y`
        pop_growth = pop_growth'
    end

    @assert length(pop_growth) == nperiods

    # Subtract log levels to get log growth rates, then take the exponential to
    # get growth rates
    100. * (exp(y./100. - y_t1./100. .+ pop_growth).^4 .- 1.)
end

"""
```
get_transform4q(transform::Function)
```
Returns the 4-quarter transformation associated with the annualizing transformation.
"""
function get_transform4q(transform::Function)

    transform4q = if transform == logtopct_annualized_percapita
        logtopct_4q_percapita
    elseif transform == logtopct_annualized
        logtopct_4q
    elseif transform == loglevelto4qpct_annualized_percapita
        loglevelto4qpct_4q_percapita
    elseif transform == quartertoannual
        quartertoannual
    elseif transform == identity
        identity
    else
        error("4q equivalent not implemented for $transform")
    end

end

"""
```
logtopct_4q(y, data, q_adj = 100)
```

Transform from log growth rates to total (not per-capita) % growth
rates (4-quarter values).
"""
function logtopct_4q(y::Array, data::Vector, q_adj = 100.)
    y = prepend_data(y, data)

    if ndims(y) == 1
        y_4q = y[1:end-3] + y[2:end-2] + y[3:end-1] + y[ 4:end]
    else
        y_4q = y[:,  1:end-3] + y[:, 2:end-2] + y[:, 3:end-1] + y[:, 4:end]
    end

    100. * (exp(y_4q/q_adj) - 1.)
end

"""
```
logtopct_4q_percapita(y, pop_growth, q_adj = 100)
```
Transform from log growth rates to 4-quarter % growth rates.

### Note

This should only be used for output, consumption, investment
and GDP deflator (inflation).

### Inputs

- `y`: The data we wish to transform to 4 quarter annualized percent change from
  1-quarter log-levels. `y` is either a vector of length `nperiods` or an
  `ndraws x `nperiods` matrix.

- `pop_growth::Vector`: The length `nperiods` vector of population growth rates.

- `q_adj`: Optional argument defaulting to 100.
"""
function logtopct_4q_percapita(y::Array, data::Vector, pop_growth::Vector, q_adj = 100)

    y = prepend_data(y, data)
    # y is now 1 x 3+H, where H is the number of forecast periods.

    pop_growth_4q = pop_growth[1:end-3] + pop_growth[2:end-2] + pop_growth[3:end-1] + pop_growth[4:end]

    # `y` is either a vector of length `nperiods` or an
    # `ndraws` x `nperiods` matrix
    if ndims(y) == 1
        y_4q = y[1:end-3] + y[2:end-2] + y[3:end-1] + y[4:end]
        nperiods = length(y_4q)
    else
        y_4q = y[:, 1:end-3] + y[:, 2:end-2] + y[:, 3:end-1] + y[:, 4:end]
        nperiods = size(y_4q, 2)

        # Transpose `pop_growth` to a 1 x `nperiods` row vector so it can be
        # broadcasted to match the dimensions of `y`
        pop_growth_4q = pop_growth_4q'
    end

    @assert length(pop_growth_4q) == nperiods

    100. * (exp(y_4q/q_adj .+ pop_growth_4q) - 1.)
end


"""
```
loglevelto4qpct_4q_percapita(y, y0, pop_growth, q_adj)
```

Transform from log level to 4-quarter percent change, adjusting for
population growth.

### Note

This is usually applied to labor supply (hours worked per hour), and
probably shouldn't be used for any other observables.

### Inputs

- `y`: The data we wish to transform to 4 quarter annualized percent change from
  1-quarter log-levels. `y` is either a vector of length `nperiods` or an
  `ndraws x `nperiods` matrix.

- `y0`: The last data point in the history (of state or observable)
  corresponding to the `y` variable.  This is required to compute a percentage
  change for the first period.

- `pop_growth::Vector`: The length `nperiods` vector of population growth rates.
"""
function loglevelto4qpct_4q_percapita{T<:AbstractFloat}(y::Array, data::Vector, pop_growth::Vector, q_adj::T = 100.)

    # `y_t4` is an array of the same size as `y`, representing the t-4
    # period observations for each draw
    nperiods = length(y)
    y_t4 = if ndims(y) == 1
        prepend_data(y[1:nperiods-4], data)
    else
        prepend_data(y[:, 1:nperiods-4], data)
    end
    y_4q = y - y_t4

    # Transpose `pop_growth` to a 1 x `nperiods` row vector so it can be
    # broadcasted to match the dimensions of `y`
    pop_growth = pop_growth'
    pop_growth_4q = pop_growth[1:end-3] + pop_growth[2:end-2] + pop_growth[3:end-1] + pop_growth[4:end]

    # Subtract log levels to get log growth rates, then exponentiate to get growth rates
    100. * (exp(y_4q./q_adj .+ pop_growth_4q) .- 1.)
end

"""
```
prepend_data(y::Array, data::Vector)
```

Prepends data necessary for running 4q transformations.

### Inputs:

- `y`: `ndraws x t` array representing a timeseries for variable `y`
- `data`: vector representing a timeseries to prepend to `y`
"""
function prepend_data(y::Array, data::Vector)
    if ndims(y) == 1
        y_extended = vcat(data, y)
    else
        ndraws = size(y, 1)
        datas  = repeat(data, ndraws, 1)
        y_extended = hcat(datas, y)
    end

    return y_extended
end
