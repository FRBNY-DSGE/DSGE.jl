function get_forecast_infiles(m, input_type)
    infiles = Dict{Symbol, AbstractString}

    if input_type == :mode
        infiles[:params] = rawpath(m,"estimate","paramsmode.h5")
    elseif input_type == :mean
        infiles[:params] = workpath(m,"estimate","paramsmean.h5")
    elseif input_type == :full
        infiles[:params] = rawpath(m,"estimate","mhsave.h5")
        infiles[:TTT]    = rawpath(m,"estimate","mhsave.h5")
        infiles[:RRR]    = rawpath(m,"estimate","mhsave.h5")
        infiles[:CCC]    = rawpath(m,"estimate","mhsave.h5")
        infiles[:zend]   = rawpath(m,"estimate","mhsave.h5")
    end

    return infiles
end
