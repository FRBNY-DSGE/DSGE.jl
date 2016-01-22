"""
`load_data(m::AbstractModel; start_date="1959-03-31", end_date=last_quarter_end())`

Checks in `inpath(m)` for vintaged datasets corresponding to the ones
in `keys(m.data_series)`. Loads the appriopriate data series
(specified in `m.data_series[source]`) for each data source, and
returns a consolidated DataFrame (merged on date).

# Arguments
- `m::AbstractModel`: the model object
- `start_date`: starting date.
- `end_date`: ending date. 

# Notes:
The name of the input data file must be the same as the source string
in `m.data_series`, and those files must be located in .txt files in
`inpath(m, "data")`.
"""
function load_data(m::AbstractModel; start_date="1959-03-31", end_date=last_quarter_end())

    # Date formatting
    start_date = Date(start_date, "y-m-d")

    if !isa(end_date, Date) end_date = Date(end_date, "y-m-d") end

    # Load FRED data, set ois series to load
    data = load_fred_data(m; start_date=firstdayofquarter(start_date), end_date=end_date)
    set_ois_series!(m)
    
    # For each additional source, search for the file with the proper name. Open
    # it, read it in, and merge it with fred_series
    
    vint = data_vintage(m)
    
    for source in keys(m.data_series)
        
        if isequal(source, :fred)
            continue
        end
        
        infile = inpath(m, "data", "$(string(source))_$vint.txt")
        mnemonics = m.data_series[source]
        
        if isfile(infile)

            # Read in dataset and check that the file contains data for the proper dates
            addl_data = readtable(infile, separator = '\t')

            # Convert dates from strings to quarter-end dates for date arithmetic
            #addl_data[:date] = stringstodates(addl_data[:date])
            #map!(x->lastdayofquarter(x), addl_data[:date], addl_data[:date])
            format_dates!(:date, addl_data)
            
            if !in(lastdayofquarter(start_date), addl_data[:date]) || 
                !in(lastdayofquarter(end_date), addl_data[:date])
                warn("$infile does not contain the entire date range specified...you may want to update your data.")
                
            end
            
            # make sure each mnemonic that was specified is present
            for series in mnemonics
                if !in(series, names(addl_data))
                    error("$(string(series)) is missing from $infile.")
                end
            end

            # extract just the columns and rows of the dataset we want, and merge them with data
            cols = [:date; mnemonics]
            rows = start_date .<= addl_data[:date] .<= end_date
            
            addl_data = addl_data[rows, cols]
            data = join(data, addl_data, on=:date, kind=:outer)
            
        elseif isempty(m.data_series[source])
            warn("No series were specificed from $(string(source))")
            continue
        else 
            error("$infile was not found." )
        end
    end

    # turn NAs into NaNs
    na2nan!(data)
    
    return sort!(data, cols = :date)
end

"""
`set_ois_series!(m::AbstractModel)`

    Sets the series to load from the file `ois_vint.txt` based on `n_anticipated_shocks(m)`.  
"""
function set_ois_series!(m::AbstractModel)
    nant = n_anticipated_shocks(m)
    if nant > 0
        m.data_series[:ois] = [symbol("ant$i") for i in 1:nant]
    end

end


"""
Returns the first day of the previous quarter
"""
function last_quarter_end()
    
    cqs = Dates.firstdayofquarter(now())     # current quarter start
    lqe = cqs - Dates.Day(1)                       # last quarter end

    Date(lqe)
end

"""
Returns a DataArray of quarter start dates between `start_date` and `end_date`.
"""
function get_quarter_ends(start_date,end_date)
    dr = start_date:end_date
    
    dates = recur(dr) do x
        lastdayofquarter(x) == x 
    end
end

"""
Converts a collection of strings in "y-m-d" format to Dates.
"""
function stringstodates(stringarray)
    n = length(stringarray)
    dates = Vector{Date}(n)
    for i = 1:n
        dates[i] = Date(stringarray[i], "y-m-d")
    end

    dates
end

function format_dates!(col, df)
    df[col] = stringstodates(df[col])
    map!(x->lastdayofquarter(x), df[col], df[col])
end

"""
```
na2nan!(df::DataFrame)
```

Convert all NAs in a DataFrame to NaNs.
"""
function na2nan!(df::DataFrame)
    for col in names(df)
        df[isna(df[col]), col] = NaN
    end
end

"""
```
na2nan!(df::DataFrame)
```

Convert all NAs in a DataFrame to NaNs.
"""
function na2nan!(v::DataArray)
    for i = 1:length(v)
        v[i] = isna(v[i]) ?  NaN : v[i]
    end
end
