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
`nominal_to_real(col, df; deflator_mnemonic = :GDPCTPI)`

Converts nominal to real values using the specified deflator.

## Arguments
- `col`: Symbol indicating which column of `df` to transform
- `df`: DataFrame containining series for proper population measure and `col`

## Keyword arguments

- `deflator_mnemonic`: indicates which deflator to use to calculate real values. Default
  value is the FRED GDP Deflator mnemonic.
"""
function nominal_to_real(col::Symbol, df::DataFrame; deflator_mnemonic::Symbol = :GDPCTPI)
    return df[col] ./ df[deflator_mnemonic]
end


"""
```
percapita(m, col, df)
percapita(col, df, population_mnemonic)
```

Converts data column `col` of DataFrame `df` to a per-capita value.

The first method checks `hpfilter_population(m)`. If true, then it divides by
the filtered population series. Otherwise it divides by the result of
`parse_population_mnemonic(m)[1]`.

## Arguments

- `col`: `Symbol` indicating which column of data to transform
- `df`: `DataFrame` containining series for proper population measure and `col`
- `population_mnemonic`: a mnemonic found in `df` for some population measure
"""
function percapita(m::AbstractModel, col::Symbol, df::DataFrame)
    if hpfilter_population(m)
        population_mnemonic = Nullable(:filtered_population)
    else
        population_mnemonic = parse_population_mnemonic(m)[1]
        if isnull(population_mnemonic)
            error("No population mnemonic provided")
        end
    end
    percapita(col, df, get(population_mnemonic))
end

function percapita(col::Symbol, df::DataFrame, population_mnemonic::Symbol)
    df[col] ./ df[population_mnemonic]
end

"""
```
series_lag_n = lag(series, n)
```
Returns a particular data series lagged by n periods
"""
function lag(series::AbstractArray, n::Int)
    return vcat(NaN*ones(n), series[1:end-n])
end

"""
```
yt, yf = hpfilter(y, λ)
```

Applies the Hodrick-Prescott filter (\"H-P filter\"). The smoothing parameter `λ` is applied
to the columns of `y`, returning the trend component `yt` and the cyclical component `yf`.
For quarterly data, one can use λ=1600.

Consecutive missing values at the beginning or end of the time series are excluded from the
filtering. If there are missing values within the series, the filtered values are all NaN.

See also:
```
Hodrick, Robert; Prescott, Edward C. (1997). \"Postwar U.S. Business Cycles: An Empirical
Investigation\". Journal of Money, Credit, and Banking 29 (1): 1–16.
```
"""
function hpfilter(y::AbstractVector, λ::Real)
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

function hpfilter_(y::AbstractVector, λ::Real)
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
difflog(x::AbstractVector)
```
"""
function difflog(x::Vector)
    [NaN; log.(x[2:end]) - log.(x[1:end-1])]
end


"""
```
difflog(x::DataArray{AbstractFloat})
```
"""
function difflog(x::DataArray)
    DSGE.na2nan!(x)
    return difflog(convert(Vector, x))
end


"""
```
oneqtrpctchange(y)
```

Calculates the quarter-to-quarter percentage change of a series.
"""
function oneqtrpctchange(y::AbstractVector)
    100 * difflog(y)
end


## REVERSE TRANSFORMS

"""
```
loggrowthtopct(y)
```

Transform from annualized quarter-over-quarter log growth rates to annualized
quarter-over-quarter percent change.

### Note

This should only be used in Model 510, which has the core PCE inflation
observable in annualized log growth rates.
"""
function loggrowthtopct(y::AbstractArray)
    100. * (exp.(y/100.) - 1.)
end

"""
```
loggrowthtopct_percapita(y, pop_growth)
```

Transform from annualized quarter-over-quarter log per-capita growth rates to
annualized quarter-over-quarter aggregate percent change.

### Note

This should only be used in Model 510, which has the output growth observable in
annualized log per-capita growth rates.

### Inputs

- `y`: the data we wish to transform to annualized percent change from
  annualized log growth rates. `y` is either a vector of length `nperiods` or an
  `ndraws x `nperiods` matrix.

- `pop_growth::Vector`: the length `nperiods` vector of log population growth
  rates.
"""
function loggrowthtopct_percapita(y::AbstractArray, pop_growth::AbstractVector)
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

    @assert length(pop_growth) == nperiods "Length of pop_growth ($(length(pop_growth))) must equal number of periods of y ($nperiods)"

    100. * ((exp.(y/100.) .* exp.(pop_growth).^4) - 1.)
end

"""
```
loggrowthtopct_annualized(y)
```

Transform from log growth rates to annualized quarter-over-quarter percent change.
"""
function loggrowthtopct_annualized(y::AbstractArray)
    100. * (exp.(y/100.).^4 - 1.)
end

"""
```
loggrowthtopct_annualized_percapita(y, pop_growth)
```

Transform from log per-capita growth rates to annualized aggregate (not
per-capita) quarter-over-quarter percent change.

### Note

This should only be used for output, consumption, investment
and GDP deflator (inflation).

### Inputs

- `y`: the data we wish to transform to annualized percent change from
  quarter-over-quarter log growth rates. `y` is either a vector of length
  `nperiods` or an `ndraws x `nperiods` matrix.

- `pop_growth::Vector`: the length `nperiods` vector of log population growth
  rates.
"""
function loggrowthtopct_annualized_percapita(y::AbstractArray, pop_growth::AbstractVector)
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

    @assert length(pop_growth) == nperiods "Length of pop_growth ($(length(pop_growth))) must equal number of periods of y ($nperiods)"

    100. * (exp.(y/100. .+ pop_growth).^4 - 1.)
end

"""
```
logleveltopct_annualized(y, y0 = NaN)
```

Transform from log levels to annualized quarter-over-quarter percent change.

### Inputs

- `y`: the data we wish to transform to annualized quarter-over-quarter percent
  change from log levels. `y` is either a vector of length `nperiods` or an
  `ndraws x `nperiods` matrix.

- `y0`: the last data point in the history (of state or observable)
  corresponding to the `y` variable. This is required to compute a percent
  change for the first period.
"""
function logleveltopct_annualized(y::AbstractArray, y0::Real = NaN)
    # `y_t1` is an array of the same size as `y`, representing the previous
    # period observations for each draw
    if ndims(y) == 1
        y_t1 = vcat([y0], y[1:end-1])
    else
        ndraws = size(y, 1)
        y0s  = fill(y0, ndraws, 1)
        y_t1 = hcat(y0s, y[:, 1:end-1])
    end

    # Subtract log levels to get log growth rates, then take the exponential to
    # get growth rates
    100. * (exp.(y./100. - y_t1./100.).^4 .- 1.)
end

"""
```
logleveltopct_annualized_percapita(y, pop_growth, y0 = NaN)
```

Transform from per-capita log levels to annualized aggregate (not per-capita)
quarter-over-quarter percent change.

### Note

This is usually applied to labor supply (hours worked per hour), and
probably shouldn't be used for any other observables.

### Inputs

- `y`: the data we wish to transform to annualized aggregate
  quarter-over-quarter percent change from per-capita log levels. `y` is either
  a vector of length `nperiods` or an `ndraws x `nperiods` matrix.

- `pop_growth::Vector`: the length `nperiods` vector of log population growth
  rates.

- `y0`: The last data point in the history (of state or observable)
  corresponding to the `y` variable. This is required to compute a percent
  change for the first period.
"""
function logleveltopct_annualized_percapita(y::AbstractArray, pop_growth::AbstractVector, y0::Real = NaN)
    # `y_t1` is an array of the same size as `y`, representing the previous
    # period observations for each draw
    if ndims(y) == 1
        nperiods = length(y)
        y_t1 = vcat([y0], y[1:end-1])
    else
        (ndraws, nperiods) = size(y)
        y0s  = fill(y0, ndraws, 1)
        y_t1 = hcat(y0s, y[:, 1:end-1])

        # Transpose `pop_growth` to a 1 x `nperiods` row vector so it can be
        # broadcasted to match the dimensions of `y`
        pop_growth = pop_growth'
    end

    @assert length(pop_growth) == nperiods "Length of pop_growth ($(length(pop_growth))) must equal number of periods of y ($nperiods)"

    # Subtract log levels to get log growth rates, then take the exponential to
    # get growth rates
    100. * (exp.(y./100. - y_t1./100. .+ pop_growth).^4 .- 1.)
end

"""
```
get_nopop_transform(transform::Function)
```
Returns the corresponding transformation which doesn't add back population
growth. Used for shock decompositions, deterministic trends, and IRFs, which are
given in deviations.
"""
function get_nopop_transform(transform::Function)
    transform4q = if transform == loggrowthtopct_annualized_percapita
        loggrowthtopct_annualized
    elseif transform == logleveltopct_annualized_percapita
        logleveltopct_annualized
    else
        transform
    end
end

"""
```
get_irf_transform(transform::Function)
```
Returns the IRF-specific transformation, which doesn't add back population
growth (since IRFs are given in deviations).
"""
function get_irf_transform(transform::Function)
    transform4q = if transform == loggrowthtopct_annualized_percapita
        loggrowthtopct_annualized
    elseif transform == logleveltopct_annualized_percapita
        logleveltopct_annualized
    else
        transform
    end
end

"""
```
get_transform4q(transform::Function)
```
Returns the 4-quarter transformation associated with the annualizing transformation.
"""
function get_transform4q(transform::Function)
    if transform == loggrowthtopct_annualized_percapita
        loggrowthtopct_4q_percapita
    elseif transform == loggrowthtopct_annualized
        loggrowthtopct_4q
    elseif transform == logleveltopct_annualized_percapita
        logleveltopct_4q_percapita
    elseif transform == logleveltopct_annualized
        logleveltopct_4q
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
loggrowthtopct_4q(y, data = fill(NaN, 3))
```

Transform from log growth rates to 4-quarter percent change.

### Inputs

- `y`: the data we wish to transform to aggregate 4-quarter percent change from
  log per-capita growth rates. `y` is either a vector of length `nperiods` or an
  `ndraws x `nperiods` matrix.

- `data`: if `y = [y_t, y_{t+1}, ..., y_{t+nperiods-1}]`, then
  `data = [y_{t-3}, y_{t-2}, y_{t-1}]`. This is necessary to compute
  4-quarter percent changes for the first three periods.
"""
function loggrowthtopct_4q(y::AbstractArray, data::AbstractVector = fill(NaN, 3))
    @assert length(data) == 3 "Length of data ($(length(data))) must be 3"

    # Prepend previous three periods to `y`
    y = prepend_data(y, data)

    # `y` is either a vector of length `nperiods+3` or an
    # `ndraws` x `nperiods+3` matrix
    if ndims(y) == 1
        y_4q = y[1:end-3] + y[2:end-2] + y[3:end-1] + y[4:end]
    else
        y_4q = y[:,  1:end-3] + y[:, 2:end-2] + y[:, 3:end-1] + y[:, 4:end]
    end

    100. * (exp.(y_4q/100.) - 1.)
end

"""
```
loggrowthtopct_4q_percapita(y, pop_growth, data = fill(NaN, 3))
```
Transform from log per-capita growth rates to aggregate 4-quarter percent
change.

### Note

This should only be used for output, consumption, investment, and GDP deflator
(inflation).

### Inputs

- `y`: the data we wish to transform to aggregate 4-quarter percent change from
  log per-capita growth rates. `y` is either a vector of length `nperiods` or an
  `ndraws x `nperiods` matrix.

- `pop_growth::Vector`: the length `nperiods` vector of log population growth
  rates.

- `data`: if `y = [y_t, y_{t+1}, ..., y_{t+nperiods-1}]`, then
  `data = [y_{t-3}, y_{t-2}, y_{t-1}]`. This is necessary to compute
  4-quarter percent changes for the first three periods.
"""
function loggrowthtopct_4q_percapita(y::AbstractArray, pop_growth::AbstractVector,
                                     data::Vector = fill(NaN, 3))
    @assert length(data) == 3 "Length of data ($(length(data))) must be 3"

    # Four-quarter population growth
    pop_growth_4q = pop_growth[1:end-3] + pop_growth[2:end-2] + pop_growth[3:end-1] + pop_growth[4:end]

    # Prepend previous three periods to `y`
    y = prepend_data(y, data)

    # `y` is either a vector of length `nperiods+3` or an
    # `ndraws` x `nperiods+3` matrix
    if ndims(y) == 1
        y_4q = y[1:end-3] + y[2:end-2] + y[3:end-1] + y[4:end]
        nperiods = length(y_4q)
    else
        y_4q = y[:, 1:end-3] + y[:, 2:end-2] + y[:, 3:end-1] + y[:, 4:end]
        nperiods = size(y_4q, 2)

        # Transpose `pop_growth` to a 1 x `nperiods` row vector so it can be
        # broadcasted to match the dimensions of `y_4q`
        pop_growth_4q = pop_growth_4q'
    end

    @assert length(pop_growth_4q) == nperiods "Length of pop_growth ($(length(pop_growth))) must equal number of periods of y ($nperiods)"

    100. * (exp.(y_4q/100. .+ pop_growth_4q) - 1.)
end

"""
```
logleveltopct_4q(y, data = fill(NaN, 4))
```

Transform from log levels to 4-quarter percent change.

### Inputs

- `y`: the data we wish to transform to 4-quarter percent change from log
  levels. `y` is either a vector of length `nperiods` or an `ndraws x `nperiods`
  matrix.

- `data`: if `y = [y_t, y_{t+1}, ..., y_{t+nperiods-1}]`, then
  `data = [y_{t-4}, y_{t-3}, y_{t-2}, y_{t-1}]`. This is necessary to compute
  4-quarter percent changes for the first three periods.
"""
function logleveltopct_4q(y::AbstractArray, data::AbstractVector = fill(NaN, 4))
    @assert length(data) == 4 "Length of data ($(length(data))) must be 4"

    # `y_t4` is an array of the same size as `y`, representing the t-4
    # period observations for each t
    y_t4 = if ndims(y) == 1
        nperiods = length(y)
        prepend_data(y[1:nperiods-4], data)
    else
        nperiods = size(y, 2)
        prepend_data(y[:, 1:nperiods-4], data)
    end
    y_4q = y - y_t4

    # Subtract log levels to get log growth rates, then exponentiate to get
    # growth rates
    100. * (exp.(y_4q./100.) .- 1.)
end

"""
```
logleveltopct_4q_percapita(y, pop_growth, data = fill(NaN, 4))
```

Transform from per-capita log levels to 4-quarter aggregate percent change.

### Note

This is usually applied to labor supply (hours worked), and probably shouldn't
be used for any other observables.

### Inputs

- `y`: the data we wish to transform to 4-quarter aggregate percent change from
  per-capita log levels. `y` is either a vector of length `nperiods` or an
  `ndraws x `nperiods` matrix.

- `pop_growth::Vector`: the length `nperiods` vector of log population growth
  rates.

- `data`: if `y = [y_t, y_{t+1}, ..., y_{t+nperiods-1}]`, then
  `data = [y_{t-4}, y_{t-3}, y_{t-2}, y_{t-1}]`. This is necessary to compute
  4-quarter percent changes for the first three periods.
"""
function logleveltopct_4q_percapita(y::AbstractArray, pop_growth::AbstractVector,
                                    data::AbstractVector = fill(NaN, 4))
    @assert length(data) == 4 "Length of data ($(length(data))) must be 4"

    # Four-quarter population growth
    pop_growth_4q = pop_growth[1:end-3] + pop_growth[2:end-2] + pop_growth[3:end-1] + pop_growth[4:end]

    # `y_t4` is an array of the same size as `y`, representing the t-4
    # period observations for each t
    if ndims(y) == 1
        nperiods = length(y)
        y_t4 = prepend_data(y[1:nperiods-4], data)
    else
        # Transpose `pop_growth` to a 1 x `nperiods` row vector so it can be
        # broadcasted to match the dimensions of `y`
        pop_growth_4q = pop_growth_4q'

        nperiods = size(y, 2)
        y_t4 = prepend_data(y[:, 1:nperiods-4], data)
    end

    @assert length(pop_growth_4q) == nperiods "Length of pop_growth ($(length(pop_growth))) must equal number of periods of y ($nperiods)"

    y_4q = y - y_t4

    # Subtract log levels to get log growth rates, then exponentiate to get growth rates
    100. * (exp.(y_4q./100. .+ pop_growth_4q) .- 1.)
end

"""
```
prepend_data(y, data)
```

Prepends data necessary for running 4q transformations.

### Inputs:

- `y`: `ndraws x t` array representing a timeseries for variable `y`
- `data`: vector representing a timeseries to prepend to `y`
"""
function prepend_data(y::AbstractArray, data::AbstractVector)
    if ndims(y) == 1
        y_extended = vcat(data, y)
    else
        ndraws = size(y, 1)
        datas  = repmat(data', ndraws, 1)
        y_extended = hcat(datas, y)
    end

    return y_extended
end

"""
```
get_scenario_transform(transform::Function)
```
Given a transformation used for usual forecasting, return the transformation
used for *scenarios*, which are forecasted in deviations from baseline.

The 1Q deviation from baseline should really be calculated by 1Q transforming
the forecasts (in levels) under the baseline (call this `y_b`) and alternative
scenario (`y_s`), then subtracting baseline from alternative scenario (since
most of our 1Q transformations are nonlinear). Let `y_d = y_s - y_b`. Then, for
example, the most correct `loggrowthtopct_annualized` transformation is:

```
y_b_1q = 100*(exp(y_b/100)^4 - 1)
y_s_1q = 100*(exp(y_s/100)^4 - 1)
y_d_1q = y_b_1q - y_s_1q
```

Instead, we approximate this by transforming the deviation directly:

```
y_d_1q ≈ 4*(y_b - y_s)
```
"""
function get_scenario_transform(transform::Function)
    if transform in [loggrowthtopct_annualized_percapita, loggrowthtopct_annualized,
                     quartertoannual]
        quartertoannual
    elseif transform in [logleveltopct_annualized_percapita, logleveltopct_annualized]
        logleveltopct_annualized_approx
    elseif transform == identity
        identity
    else
        error("Scenario equivalent not implemented for $transform")
    end
end

function get_scenario_transform4q(transform::Function)
    if transform in [loggrowthtopct_annualized_percapita, loggrowthtopct_annualized]
        loggrowthtopct_4q_approx
    elseif transform in [logleveltopct_annualized_percapita, logleveltopct_annualized]
        logleveltopct_4q_approx
    elseif transform == quartertoannual
        quartertoannual
    elseif transform == identity
        identity
    else
        error("Scenario 4q equivalent not implemented for $transform")
    end
end

"""
```
logleveltopct_annualized_approx(y, y0 = NaN)
```

Transform from log levels to *approximate* annualized quarter-over-quarter
percent change.

**This method should only be used to transform scenarios forecasts, which are in
  deviations from baseline.**

### Inputs

- `y`: the data we wish to transform to annualized quarter-over-quarter percent
  change from log levels. `y` is either a vector of length `nperiods` or an
  `ndraws x `nperiods` matrix.

- `y0`: the last data point in the history (of state or observable)
  corresponding to the `y` variable. This is required to compute a percent
  change for the first period.
"""
function logleveltopct_annualized_approx(y::AbstractArray, y0::Real = NaN)
    # `y_t1` is an array of the same size as `y`, representing the previous
    # period observations for each draw
    if ndims(y) == 1
        y_t1 = vcat([y0], y[1:end-1])
    else
        ndraws = size(y, 1)
        y0s  = fill(y0, ndraws, 1)
        y_t1 = hcat(y0s, y[:, 1:end-1])
    end

    # Subtract log levels to get log growth rates, then multiply by 4 to
    # approximate annualizing
    4(y - y_t1)
end

"""
```
loggrowthtopct_4q_approx(y, data = fill(NaN, 3))
```

Transform from log growth rates to *approximate* 4-quarter percent change.

**This method should only be used to transform scenarios forecasts, which are in
  deviations from baseline.**

### Inputs

- `y`: the data we wish to transform to aggregate 4-quarter percent change from
  log per-capita growth rates. `y` is either a vector of length `nperiods` or an
  `ndraws x `nperiods` matrix.

- `data`: if `y = [y_t, y_{t+1}, ..., y_{t+nperiods-1}]`, then
  `data = [y_{t-3}, y_{t-2}, y_{t-1}]`. This is necessary to compute
  4-quarter percent changes for the first three periods.
"""
function loggrowthtopct_4q_approx(y::AbstractArray, data::AbstractVector = fill(NaN, 3))
    @assert length(data) == 3 "Length of data ($(length(data))) must be 3"

    # Prepend previous three periods to `y`
    y = prepend_data(y, data)

    # `y` is either a vector of length `nperiods+3` or an
    # `ndraws` x `nperiods+3` matrix
    if ndims(y) == 1
        y[1:end-3] + y[2:end-2] + y[3:end-1] + y[4:end]
    else
        y[:,  1:end-3] + y[:, 2:end-2] + y[:, 3:end-1] + y[:, 4:end]
    end
end

"""
```
logleveltopct_4q_approx(y, data = fill(NaN, 4))
```

Transform from log levels to *approximate* 4-quarter percent change.

**This method should only be used to transform scenarios forecasts, which are in
  deviations from baseline.**

### Inputs

- `y`: the data we wish to transform to 4-quarter percent change from log
  levels. `y` is either a vector of length `nperiods` or an `ndraws x `nperiods`
  matrix.

- `data`: if `y = [y_t, y_{t+1}, ..., y_{t+nperiods-1}]`, then
  `data = [y_{t-4}, y_{t-3}, y_{t-2}, y_{t-1}]`. This is necessary to compute
  4-quarter percent changes for the first three periods.
"""
function logleveltopct_4q_approx(y::AbstractArray, data::AbstractVector = fill(NaN, 4))
    @assert length(data) == 4 "Length of data ($(length(data))) must be 4"

    # `y_t4` is an array of the same size as `y`, representing the t-4
    # period observations for each t
    y_t4 = if ndims(y) == 1
        nperiods = length(y)
        prepend_data(y[1:nperiods-4], data)
    else
        nperiods = size(y, 2)
        prepend_data(y[:, 1:nperiods-4], data)
    end

    y - y_t4
end
