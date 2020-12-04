"""
```
load_fred_data(m::AbstractDSGEModel; start_date="1959-03-31", end_date=prev_quarter())
```

Checks in `inpath(m, raw)` for a FRED dataset corresponding to `data_vintage(m)`.
If a FRED vintage exists on disk, any required FRED series that is contained therein will be
imported. All missing series will be downloaded directly from FRED using the *FredData*
package. The full dataset is written to the appropriate data vintage file and returned.

# Arguments
- `m::AbstractDSGEModel`: the model object
- `start_date`: starting date.
- `end_date`: ending date.

# Notes

The FRED API reports observations according to the quarter-start date. `load_fred_data`
returns data indexed by quarter-end date for compatibility with other datasets.
"""
function load_fred_data(m::AbstractDSGEModel;
                        start_date::Dates.Date = Dates.Date("1959-01-01", "y-m-d"),
                        end_date::Dates.Date   = prev_quarter(),
                        verbose::Symbol  = :low)

    data_series = parse_data_series(m)
    if haskey(data_series, :FRED)
        mnemonics = data_series[:FRED]
    else
        # Then no FRED data used, return DataFrame with dates
        return DataFrame(date = get_quarter_ends(start_date,end_date))
    end
    vint = data_vintage(m)
    dateformat = "yymmdd"

    # Set up dataset and labels
    missing_series = Vector{Symbol}()
    data = []

    datafile = inpath(m, "raw", "fred_$vint.csv")
    if isfile(datafile)

        # Read in dataset and check that the file contains data for the proper dates
        data = CSV.read(datafile, DataFrame)

        # Convert dates from strings to dates for date arithmetic
        format_dates!(:date, data)

        qstart = lastdayofquarter(start_date)
        qend = lastdayofquarter(end_date)

        if !in(qstart, data[!,:date]) || !in(qend, data[!,:date])
            println(verbose, :low, "FRED dataset on disk missing start or end date. Fetching data from FRED.")
            data = DataFrame(date = get_quarter_ends(start_date,end_date))
            missing_series = mnemonics
        else
            # If we have the right dates but the series isn't in the file, add it to
            # missing_series
            for series in mnemonics
                if !in(series, propertynames(data))
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

        # Have to do this wacky parsing to prepend the century to the data vintage
        vint_date = if parse(Int,vint[1:2]) < 59
            Year(2000) + Dates.Date(vint, dateformat)
        else
            Year(1900) + Dates.Date(vint, dateformat)
        end

        fredseries = Array{FredSeries, 1}(undef, length(missing_series))
        f = Fred()

        for (i,s) in enumerate(missing_series)
            println(verbose, :low, "Fetching FRED series $s...")
            try
                fredseries[i] = get_data(f, string(s); frequency="q",
                                                       observation_start=string(start_date),
                                                       observation_end=string(end_date),
                                                       vintage_dates=string(vint_date))
            catch err
                if :msg in fieldnames(typeof(err))
                    @warn err.msg
                else
                    show(err)
                end
                @warn "FRED series $s could not be fetched at vintage $vint."

                try
                    println(verbose, :low, "Fetching FRED series $s without vintage...")
                    fredseries[i] = get_data(f, string(s); frequency="q",
                                                           observation_start=string(start_date),
                                                           observation_end=string(end_date))
                catch err
                    if :msg in fieldnames(typeof(err))
                        @warn err.msg
                    else
                        show(err)
                    end
                    @warn "FRED series $s could not be fetched."
                    continue
                end
            end
        end

        # Extract dataframe from each series and merge on date
        has_oj = isdefined(DataFrames, :outerjoin)
        for i = 1:length(fredseries)
            if isassigned(fredseries, i)
                series = fredseries[i]
                series_id = Symbol(series.id)
                rename!(series.df, :value => series_id)
                map!(x->lastdayofquarter(x), series.df[!, :date], series.df[!, :date])
                data = has_oj ? outerjoin(data, series.df[!, [:date, series_id]], on = :date) :
                    join(data, series.df[!, [:date, series_id]], on = :date, kind = :outer)
            end
        end

        # Change the dates to be the last day of each quarter
        n_rows, n_cols = size(data)

        for i = 1:n_rows
            data[i,:date] = Dates.lastdayofquarter(data[i,:date])
        end

        if !m.testing
            CSV.write(datafile, data, missingstring = "")
            println(verbose, :low, "Updated data from FRED written to $datafile.")
        end
    end

    # Make sure to only return the series and dates that are specified for this
    # model (there may be additional series in the file)

    rows = start_date .<= data[!, :date] .<= end_date
    cols = [:date; mnemonics]
    return data[rows, cols]
end
