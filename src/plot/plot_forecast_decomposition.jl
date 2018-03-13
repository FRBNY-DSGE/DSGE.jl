function make_forecast_decomposition_mb(m_new::AbstractModel, m_old::AbstractModel,
                                        input_type::Symbol, cond_new::Symbol, cond_old::Symbol, class::Symbol)
    # Read in means
    input_file = get_decomp_mean_file(m_new, m_old, input_type, cond_new, cond_old, class)
    decomps = jldopen(input_file, "r") do file
        read(file, "decomps")
    end

    # Common metadata
    comps = [:state, :shock, :data, :param]
    dates = decomps[collect(keys(decomps))[1]][:date]
    vars  = collect(keys(DSGE.get_dict(m_new, class)))

    metadata = Dict{Symbol, Any}()
    metadata[:para]            = input_type
    metadata[:indices]         = DSGE.get_dict(m_new, class)
    metadata[:class]           = class
    metadata[:date_inds]       = DataStructures.OrderedDict(date => i for (i, date) in enumerate(dates))
    metadata[:forecast_string] = ""
    metadata[:cond_type]       = cond_new

    # Shock decomposition
    shockdec_mb = MeansBands(Dict(metadata), DataFrame(date = dates), Dict{Symbol, DataFrame}())
    shockdec_mb.metadata[:product] = :shockdec
    shockdec_mb.metadata[:shock_indices] = DataStructures.OrderedDict(comp => i for (i, comp) in enumerate(comps))
    for comp in comps
        for var in vars
            varcomp = Symbol(var, "__", comp)
            shockdec_mb.means[varcomp] = decomps[var][comp]
            shockdec_mb.bands[varcomp] = DataFrame(date = dates)
        end
    end

    # Trend
    trend_mb = MeansBands(Dict(metadata), DataFrame(date = dates), Dict{Symbol, DataFrame}())
    trend_mb.metadata[:product] = :trend
    for var in vars
        trend_mb.means[var] = 0.0
        trend_mb.bands[var] = DataFrame(date = dates)
    end

    # Deterministic trend
    dettrend_mb = MeansBands(Dict(metadata), DataFrame(date = dates), Dict{Symbol, DataFrame}())
    dettrend_mb.metadata[:product] = :dettrend
    for var in vars
        dettrend_mb.means[var] = 0.0
        dettrend_mb.bands[var] = DataFrame(date = dates)
    end

    # History
    hist_inds  = dates .<= date_mainsample_end(m_new)
    hist_dates = dates[hist_inds]
    if isempty(hist_dates)
        hist_mb = MeansBands()
    else
        hist_mb = MeansBands(Dict(metadata), DataFrame(date = hist_dates), Dict{Symbol, DataFrame}())
        hist_mb.metadata[:product]   = :hist
        hist_mb.metadata[:date_inds] = DataStructures.OrderedDict(date => i for (i, date) in enumerate(hist_dates))
        for var in vars
            hist_mb.means[var] = decomps[var][hist_inds, :total]
            hist_mb.bands[var] = DataFrame(date = dates)
        end
    end

    # Forecast
    fcast_inds  = dates .> date_mainsample_end(m_new)
    fcast_dates = dates[fcast_inds]
    if isempty(fcast_dates)
        fcast_mb = MeansBands()
    else
        fcast_mb = MeansBands(Dict(metadata), DataFrame(date = fcast_dates), Dict{Symbol, DataFrame}())
        fcast_mb.metadata[:product]   = :forecast
        fcast_mb.metadata[:date_inds] = DataStructures.OrderedDict(date => i for (i, date) in enumerate(fcast_dates))
        for var in vars
            fcast_mb.means[var] = decomps[var][fcast_inds, :total]
            fcast_mb.bands[var] = DataFrame(date = dates)
        end
    end

    return shockdec_mb, trend_mb, dettrend_mb, hist_mb, fcast_mb
end

function plot_forecast_decomposition(m_new::AbstractModel, m_old::AbstractModel,
                                     vars::Vector{Symbol}, class::Symbol,
                                     input_type::Symbol, cond_new::Symbol, cond_old::Symbol;
                                     titles::Vector{String} = String[],
                                     verbose::Symbol = :low,
                                     kwargs...)
    # Create MeansBands
    mbs = make_forecast_decomposition_mb(m_new, m_old, input_type, cond_new, cond_old, class)

    # Create shock grouping
    groups = [ShockGroup("state", [:state], :cyan),
              ShockGroup("shock", [:shock], :orange),
              ShockGroup("data",  [:data],  :green),
              ShockGroup("param", [:param], :purple)]

    # Get titles if not provided
    if isempty(titles)
        detexify_title = typeof(Plots.backend()) == Plots.GRBackend
        titles = map(var -> DSGE.describe_series(m_new, var, class, detexify = detexify_title), vars)
    end

    # Loop through variables
    plots = DataStructures.OrderedDict{Symbol, Plots.Plot}()
    for (var, title) in zip(vars, titles)
        # Call recipe
        plots[var] = DSGE.shockdec(var, mbs..., groups;
                              ylabel = DSGE.series_ylabel(m_new, var, class),
                              title = title, kwargs...)

        # Save plot
        output_file = get_decomp_filename(m_new, m_old, input_type, cond_new, cond_old, :decomp, class,
                                          pathfcn = figurespath, fileformat = DSGE.plot_extension())

        DSGE.save_plot(plots[var], output_file, verbose = verbose)
    end
end