"""
```
make_decomp_mbs(m_new, m_old, input_type, cond_new, cond_old, class, hs;
    individual_shocks = false)
```

Construct and return the `MeansBands` (for shockdec, trend, dettrend, hist, and
forecast) necessary to call the plotting function `shockdec` in
`plot_forecast_decomposition`.
"""
function make_decomp_mbs(m_new::M, m_old::M, input_type::Symbol,
                         cond_new::Symbol, cond_old::Symbol,
                         class::Symbol, hs::UnitRange{Int};
                         individual_shocks::Bool = false) where M<:AbstractModel
    # Read in means
    input_file = get_decomp_mean_file(m_new, m_old, input_type, cond_new, cond_old, class, hs)
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
    metadata[:date_inds]       = OrderedDict(date => i for (i, date) in enumerate(dates))
    metadata[:forecast_string] = ""
    metadata[:cond_type]       = cond_new

    # Shock decomposition
    shockdec_mb = MeansBands(Dict(metadata), DataFrame(date = dates), Dict{Symbol, DataFrame}())
    shockdec_mb.metadata[:product] = :shockdec
    if individual_shocks
        for var in vars
            shocks = setdiff(names(decomps[var]), vcat([:date, :total], comps))
            if var == vars[1]
                shockdec_mb.metadata[:shock_indices] = OrderedDict(shock => i for (i, shock) in enumerate(shocks))
            end
            for shock in shocks
                varshock = Symbol(var, "__", shock)
                shockdec_mb.means[varshock] = decomps[var][shock]
                shockdec_mb.bands[varshock] = DataFrame(date = dates)
            end
        end
    else
        shockdec_mb.metadata[:shock_indices] = OrderedDict(comp => i for (i, comp) in enumerate(comps))
        for var in vars
            for comp in comps
                varcomp = Symbol(var, "__", comp)
                shockdec_mb.means[varcomp] = decomps[var][comp]
                shockdec_mb.bands[varcomp] = DataFrame(date = dates)
            end
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
        hist_mb.metadata[:date_inds] = OrderedDict(date => i for (i, date) in enumerate(hist_dates))
        for var in vars
            hist_mb.means[var] = decomps[var][hist_inds, individual_shocks ? :shock : :total]
            hist_mb.bands[var] = DataFrame(date = hist_dates)
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
        fcast_mb.metadata[:date_inds] = OrderedDict(date => i for (i, date) in enumerate(fcast_dates))
        for var in vars
            fcast_mb.means[var] = decomps[var][fcast_inds, individual_shocks ? :shock : :total]
            fcast_mb.bands[var] = DataFrame(date = fcast_dates)
        end
    end

    return shockdec_mb, trend_mb, dettrend_mb, hist_mb, fcast_mb
end

"""
```
plot_forecast_decomposition(m_new, m_old, vars, class, input_type,
    cond_new, cond_old, hs; titles = [], individual_shocks = false,
    groups = shock_groupings(m_new), verbose = :low, kwargs...)
```

Plot forecast decomposition (looks like a shock decomposition). If
`individual_shocks = false`, then the black and red lines give the total
difference and the bars give the state, shock, data, and parameter
components. Otherwise, the black and red lines give the total shock component
and the bars give the individual shock contributions.

The `groups` keyword argument is only used if `individual_shocks = true`.
"""
function plot_forecast_decomposition(m_new::M, m_old::M, vars::Vector{Symbol}, class::Symbol,
                                     input_type::Symbol, cond_new::Symbol,
                                     cond_old::Symbol, hs::UnitRange{Int};
                                     titles::Vector{String} = String[],
                                     individual_shocks::Bool = false,
                                     groups::Vector{ShockGroup} = shock_groupings(m_new),
                                     verbose::Symbol = :low,
                                     kwargs...) where M<:AbstractModel
    # Create MeansBands
    mbs = make_decomp_mbs(m_new, m_old, input_type, cond_new, cond_old, class, hs,
                          individual_shocks = individual_shocks)

    # Create shock grouping
    if individual_shocks
        # Remove deterministic trend (zero by construction)
        ind_dt = findfirst(group -> group.name == "dt", groups)
        splice!(groups, ind_dt)
    else
        groups = [ShockGroup("state", [:state], colorant"#E5FCC2"), # yellow-green
                  ShockGroup("shock", [:shock], colorant"#9DE0AD"), # sea foam green
                  ShockGroup("data",  [:data],  colorant"#45ADA8"), # turquoise
                  ShockGroup("param", [:param], colorant"#547980")] # blue gray
    end

    # Get titles if not provided
    if isempty(titles)
        detexify_title = typeof(Plots.backend()) == Plots.GRBackend
        titles = map(var -> DSGE.describe_series(m_new, var, class, detexify = detexify_title), vars)
    end

    # Loop through variables
    plots = OrderedDict{Symbol, Plots.Plot}()
    for (var, title) in zip(vars, titles)
        # Call recipe
        plots[var] = DSGE.shockdec(var, mbs..., groups;
                              ylabel = DSGE.series_ylabel(m_new, var, class),
                              title = title, kwargs...)

        # Save plot
        basename = Symbol(:decomp, individual_shocks ? :shocks : :total, "_", var)
        output_file = get_decomp_filename(m_new, m_old, input_type, cond_new, cond_old,
                                          basename, Symbol(), hs,
                                          pathfcn = figurespath, fileformat = DSGE.plot_extension())

        DSGE.save_plot(plots[var], output_file, verbose = verbose)
    end
end