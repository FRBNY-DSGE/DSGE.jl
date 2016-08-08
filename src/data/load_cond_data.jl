function load_cond_data(m, cond_type)
    # Preapare file name
    cond_vintage = get_setting(m, :cond_vintage)
    cond_id = get_setting(m, :cond_id)
    cond_file_name = inpath(m, "cond", "cond_vint=$(cond_vintage)_cdid=$(cond_id).csv")

    # Read data
    cond_data = readtable(cond_file_name)
    format_dates!(:date, cond_data)

    # Make sure each observable in `m.observables` is present
    for obs in keys(m.observables)
        if !in(obs, names(cond_data) )
            error("$(string(obs)) is missing from $cond_file_name.") 
        end 
    end 

    # If we are using semi-conditional data only, then we must index out the
    # semi-conditional columns and replace all other columns with NaN.
    if cond_type == :semi
        cond_semi_names = get_setting(m, :cond_semi_names)
        cond_semi_names_nan = setdiff(names(cond_data), push!(cond_semi_names, :date))
        cond_data[:, cond_semi_names_nan] = convert(eltype(cond_data), NaN)
    end

    return cond_data
end
