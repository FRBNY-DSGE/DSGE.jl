function load_cond_data(m, cond_type; verbose::Symbol=:low)
    
    if cond_type == :none
        return DataFrame()
    end

    # Prepare file name
    cond_vintage = get_setting(m, :cond_vintage)
    cond_id = get_setting(m, :cond_id)
    file = inpath(m, "cond", "cond_vint=$(cond_vintage)_cdid=$(cond_id).csv")

    if isfile(file)
        if VERBOSITY[verbose] >= VERBOSITY[:low]
            println("Reading conditional data from $file...")
        end

        # Read data
        cond_data = readtable(file)
        format_dates!(:date, cond_data)

        # Make sure each mnemonic that was specified is present
        mnemonics = m.data_series[:conditional]
        for series in mnemonics
            if !in(series, names(cond_data))
                error("$(string(series)) is missing from $file.")
            end
        end

        return cond_data
    else
        # If series not found, throw an error
        error("Conditional data in $file not found")
    end
end
