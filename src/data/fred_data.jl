"""
`load_fred_data(m::AbstractModel; start_date="1959-01-01", end_date=last_quarter_end())`

Checks in `inpath(m)` for a FRED dataset corresponding to
`data_vintage(m)`. Necessary series (identified by series mnemonics in
`m.data_series[:fred]`) that are not found on disk are fetched from the FRED
API. The full dataset is returned and written to the appropriate data vintage file.

# Arguments
- `m::AbstractModel`: the model object
- `start_date`: starting date. 
"""
function load_fred_data(m::AbstractModel; start_date="1959-01-01", end_date=last_quarter_end())
    mnemonics = m.data_series[:fred]
    
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
    
    datafile = inpath(m, "data", "data_$vint.txt") 
    if isfile(datafile)

        # Read in dataset and check that the file contains data for the proper dates
        dataset = readtable(datafile, separator = '\t')

        if !in(string(firstdayofquarter(start_date)), dataset[:date]) || 
            !in(string(firstdayofquarter(end_date)), dataset[:date])
            missing_series = mnemonics
        else
            # If we have the right dates but the series isn't in the file, add it to missing_series
            for series in mnemonics
                if !in(symbol(series), names(dataset))
                    append!(missing_series,series)
                end
            end
        end
    else
        missing_series = mnemonics
        dataset = DataFrame(date = get_quarter_starts(start_date,end_date))
    end

    # Get the missing data series
    if missing_series != []
        
        fredseries = Array{FredSeries, 1}(length(missing_series))
        f = Fred()
        
        for (i,s) in enumerate(missing_series)
            println(s)
            try
                fredseries[i] = get_data(f, s; frequency="q", observation_start = string(start_date),
                                         observation_end=string(end_date), vintage_dates=string(vint_date))
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
                dataset = join(dataset, series.df[:date, series_id], on=:date)
            end
        end

        writetable(datafile, dataset, separator = '\t')
    end

    # Make sure to only return the series and dates that are specified for this
    # model (there may be additional series in the # file)

     cols = [:date; convert(Vector{Symbol}, [symbol(s) for s in mnemonics])]
     return dataset[string(start_date):string(end_date), cols]
end
