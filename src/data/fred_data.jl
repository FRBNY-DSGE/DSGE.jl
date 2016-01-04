"""
`load_fred_data(m::AbstractModel; start_date="1959-03-31", end_date=last_quarter_end())`

Checks in `inpath(m)` for a FRED dataset corresponding to
`data_vintage(m)`. The appropriate vintage of necessary series
(identified by series mnemonics in `m.data_series[:fred]`) that are
not found on disk are fetched via the FRED API. The full dataset is
written to the appropriate data vintage file and returned.

# Arguments
- `m::AbstractModel`: the model object
- `start_date`: starting date.
- `end_date`: ending date. 

# Notes

The FRED API reports observations according to the quarter-start
date. `load_fred_data` returns data indexed by quarter-end date for
compatibility with other datasets.
"""
function load_fred_data(m::AbstractModel;
                        start_date::Date=Date("1959-01-01", "y-m-d"),
                        end_date::Date=last_quarter_end())
    
    mnemonics = m.data_series[:fred]
    vint = data_vintage(m)
    df2 = "yymmdd"

    # Have to do this wacky parsing to prepend the century to the data vintage
    vint_date = if parse(Int,vint[1:2]) < 59
        Year(2000) + Date(vint, df2)
    else
        Year(1900) + Date(vint, df2)
    end

    # Get data vintage and check that last_quarter_end is before that.
    @assert vint_date >= end_date "Data vintage must be more recent than end date specified."
    
    n_quarters = 4*(year(end_date) - year(start_date)) - quarterofyear(start_date) + quarterofyear(end_date) + 1
    
    # Set up dataset and labels
    missing_series = Vector{Symbol}()
    data = []

    datafile = inpath(m, "data", "fred_$vint.txt")
    if isfile(datafile)
        
        # Read in dataset and check that the file contains data for the proper dates
        data = readtable(datafile, separator = '\t')

        # Convert dates from strings to dates for date arithmetic
        data[:date] = stringstodates(data[:date])
        println(names(data))
        map!(x->lastdayofquarter(x), data[:date], data[:date])
        
        qstart = lastdayofquarter(start_date)
        qend = lastdayofquarter(end_date)

        if !in(qstart, data[:date]) || !in(qend, data[:date])
            println("Start and end dates are not in data...fetching all FRED series...")
            data = DataFrame(date = get_quarter_ends(start_date,end_date))
            missing_series = mnemonics
        else
            # If we have the right dates but the series isn't in the file, add it to missing_series
            for series in mnemonics
                if !in(series, names(data))
                    push!(missing_series,series)
                end
            end
        end
    else
        missing_series = mnemonics
        data = DataFrame(date = get_quarter_ends(start_date,end_date))
    end

    # Get the missing data series
    if !isempty(missing_series)
        
        fredseries = Array{FredSeries, 1}(length(missing_series))
        f = Fred()
        
        for (i,s) in enumerate(missing_series)
            println("Fetching $(s)...")
            try
                fredseries[i] = get_data(f, string(s); frequency="q", observation_start=string(start_date),
                                         observation_end=string(end_date), vintage_dates=string(vint_date))
            catch err
                warn("$s could not be fetched.")
                continue
            end
        end
        
        # Extract dataframe from each series and merge on date
        for i = 1:length(fredseries)
            if isdefined(fredseries, i)
                series = fredseries[i]
                series_id = symbol(series.id)
                rename!(series.df, :value, series_id)
                map!(x->lastdayofquarter(x), series.df[:date], series.df[:date])
                data = join(data, series.df[:,[:date, series_id]], on=:date)
            end
        end
        
        # Change the dates to be the last day of each quarter
        n_rows, n_cols = size(data)

        for i = 1:n_rows
            data[i,:date] = Dates.lastdayofquarter(data[i,:date])
        end
        
        writetable(datafile, data, separator = '\t')
        println("Data written to $datafile.")
    end

    # Make sure to only return the series and dates that are specified for this
    # model (there may be additional series in the file)

    rows = start_date .<= data[:date] .<= end_date
    cols = [:date; mnemonics]
    return data[rows, cols]
end
