"""
`load_fred_data(m::AbstractModel; start_date="1959-01-01", end_date=last_quarter_end())`

Checks in `inpath(m)` for a FRED dataset corresponding to
`data_vintage(m)`. Necessary series (identified by series mnemonics in
`m.fred_series`) that are not found on disk are fetched from the FRED
API. The full dataset is returned and written to the appropriate data vintage file.
"""
function load_fred_data(m::AbstractModel; start_date="1959-01-01", end_date=last_quarter_end())
    
    # Date formatting
    df = "y-m-d"
    start_date = Date(start_date, df)
    if isa(end_date, AbstractString) end_date = Date(end_date, df) end

    # Get data vintage and check that last_quarter_end is before that.
    vint = data_vintage(m)
    vint_date = parse(Int,vint[1:2]) < 59 ? Year(2000) + Date(vint, "yymmdd") : Year(1900) + Date(vint, "yymmdd")
    vint_date < end_date && error("Data vintage is older than the end date specified")
    
    n_quarters = 4*(year(end_date) - year(start_date)) - quarterofyear(start_date) + quarterofyear(end_date) + 1
    
    # Set up dataset and labels
    missing_series = []
    dataset = []
    
    # Check to see if a dataset of our particular vintage exists and that we have the right number of observations
    datafile = inpath(m, "data", "data_$vint.txt") 
    if isfile(datafile)

        dataset = readtable(datafile, separator = '\t')
        @assert size(dataset)[1] == n_quarters "Number of quarters requested does not match number of quarters of data on disk."
        
        # If series isn't in the file, add it to missing_series
        for series in m.fred_series
            if !in(symbol(series), names(dataset))
                append!(missing_series,series)
            end
        end

    else
        missing_series = m.fred_series
        dataset = DataFrame(date = get_quarter_starts(start_date,end_date))
    end

    # Get the missing data series
    if missing_series != []

        fredseries = Array{FredSeries, 1}(length(missing_series))
        f = Fred()

        for (i,s) in enumerate(missing_series)
            println(s)
            try
                fredseries[i] = get_data(f, s; frequency="q", observation_start=string(start_date), observation_end=string(end_date))
            catch err
                warn("$s could not be fetched.")
                continue
            end
        end
        
        # Extract dataframes from series and merge them on date
        for i = 1:length(fredseries)
            if isdefined(fredseries, i)
                series = fredseries[i]
                series_id = symbol(series.id)
                rename!(series.df, :value, series_id)
                dataset = join(dataset, series.df[[:date, series_id]], on=:date)
            end
        end

    end


    writetable(datafile, dataset, separator = '\t')
    return dataset
end



"""
Returns the first day of the previous quarter
"""
function last_quarter_end()
    
    cqs = Dates.firstdayofquarter(now())     # current quarter start
    lqe = cqs - Day(1)                       # last quarter end

    Dates.format(lqe, "yyyy-mm-dd")
end

"""
Returns a DataArray of quarter start dates between `start_date` and `end_date`.
"""
function get_quarter_starts(start_date,end_date)
    dr = start_date:end_date
    
    dates = recur(dr) do x
        firstdayofquarter(x) == x 
    end
end
