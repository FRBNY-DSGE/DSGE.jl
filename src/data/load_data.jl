function load_data(m::AbstractModel; start_date="1959-01-01", end_date=last_quarter_end())
    data = load_fred_data(m, start_date, end_date)

    # For each dataset, search for the file with the proper name. Open
    # it, read it in, and merge it with fred_series

    vint = data_vintage(m)
    for source in keys(m.data_series)
        infile = inpath(m, "data", "$(string(source))_$vint")
        mnemonics = m.data_series[source]
        dataset = []
        
        if isfile(infile)

            # Read in dataset and check that the file contains data for the proper dates
            addl_data = readtable(datafile, separator = '\t')
            
            if !in(string(firstdayofquarter(start_date)), dataset[:date]) || 
                !in(string(firstdayofquarter(end_date)), dataset[:date])
                error("$infile does not contain the date range specified. You may need to update your data.")
            end

            # make sure each mnemonic that was specified is present
            for series in mnemonics
                if !in(symbol(series), names(dataset))
                    error("$(series) is missing from $infile.")
                end
            end

            # extract just the columns and rows of the dataset we want, and merge them with data
            cols = [:date; convert(Vector{Symbol}, [symbol(s) for s in mnemonics])]]
            addl_data = addl_data[start_date:string(end_date), cols]            
            data = join(data, addl_data, on=:date)
            
        else
            error("$infile was not found." )
        end

    end

    return data
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

