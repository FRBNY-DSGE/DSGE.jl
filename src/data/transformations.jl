function Base.call(x::Function, df::DataFrame)
    x(df)
end


"""
`function annualtoquarter(v)`

Convert from annual to quarter by dividing by 4
"""
function annualtoquarter(v)
    v / 4
end

"""
`nominal_to_real(col, df; deflator=:GDP)`

Converts nominal to real values using the specified
deflator.

## Parameters
- `col`: symbol indicating which column of `df` to transform
- `df`: DataFrame containining series for proper population measure and `col`

## Keyword arguments

- `deflator`: indicates which deflator to use to calculate real values.
   - `:GDP`: GDP deflator
       - FRED mnemonic: :GDPCTPI
   - `:PCE`: PCE deflator
       - FRED mnemonic: :PCEPILFE
"""
function nominal_to_real(col, df; deflator=:GDP)
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

Convert column `col` in `df` from a nominal to a real, per-person-employed value
"""
function nominal_to_realpercapita(col, df; population_measure=:adult_population, deflator=:GDP, scale = 1)

    # create a temporary column to hold per-capita values
    df[:temp] = percapita(col, df, population_measure=population_measure) 

    # convert to real values and scale 
    realpercapita = scale * nominal_to_real(:temp, df, deflator=deflator)

    # delete temporary column
    delete!(df, :temp)
    
    return realpercapita
end


"""
`percapita(col, df; population_measure=:adult_population)`

Converts data column `col` of DataFrame `df` to a per-capita value.

## Parameters
- `col`: symbol indicating which column of data to transform
- `df`: DataFrame containining series for proper population measure and `col`

## Keyword arguments
- `population_measure`: currently implemented for the following 3 values:
   - `:adult_population`
      - FRED mnemonic: CNP16OV
      - description: civilian noninstitutional population age 16+
      - source: Bureau of Labor Statistics
   - `:nonfarm_employed`
      - FRED mnemonic: :PRS85006013
      - description: nonfarm business sector employment
      - source: Bureau of Labor Statistics
   - `:civilian_employed`
      - FRED mnemonic: :CE16OV
      - description: civilian employment
      - source: Bureau of Labor Statistics
"""
function percapita(col, df; population_measure=:adult_population)
    
    if population_measure == :adult_population
        df[col] ./ df[:CNP16OV]
    elseif population_measure == :nonfarm_employed
        df[col] ./ df[:PRS85006013]
    elseif population_measure == :civilian_employed
        df[col] ./ df[:CE16OV]
    else
        error("percapita isn't defined for this population measure.")
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
        a[i,i]   = 6λ+1;
        a[i,i+1] = -4λ;
        a[i,i+2] = λ;
        a[i,i-2] = λ;
        a[i,i-1] = -4λ;
    end

    a[2,2] = 1+5λ;
    a[2,3] = -4λ;
    a[2,4] = λ;
    a[2,1] = -2λ;
    a[1,1] = 1+λ;
    a[1,2] = -2λ;
    a[1,3] = λ;
    
    a[n-1,n-1] = 1+5λ;
    a[n-1,n-2] = -4λ;
    a[n-1,n-3] = λ;
    a[n-1,n]   = -2λ;
    a[n,n]     = 1+λ;
    a[n,n-1]   = -2λ;
    a[n,n-2]   = λ;

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
