function quarter_date_to_number(date::Date)
    y = Dates.year(date)
    q = Dates.quarterofyear(date)
    return y + (q-1)*0.25
end

function quarter_number_to_date(datenum::Real)
    if datenum % 0.25 != 0
        throw(DomainError())
    end

    y = convert(Int, floor(datenum))
    q = convert(Int, (datenum % 1) / 0.25) + 1
    return quartertodate("$y-Q$q")
end

function plot_extension()
    be = typeof(Plots.backend())
    if be == Plots.GRBackend
        :pdf
    elseif be in [Plots.PlotlyBackend, Plots.PlotlyJSBackend]
        :html
    else
        :pdf
    end
end

function describe_series(m::AbstractModel, var::Symbol, class::Symbol;
                         detexify::Bool = false)
    res = if class in [:obs, :pseudo]
        dict = if class == :obs
            m.observable_mappings
        elseif class == :pseudo
            m.pseudo_observable_mappings
        end
        dict[var].name
    elseif class == :states
        string(var)
    elseif class in [:shocks, :stdshocks]
        replace(string(var), r"_sh$", "")
    else
        error("Invalid class: " * string(class))
    end

    detexify ? DSGE.detexify(res) : res
end

function series_ylabel(m::AbstractModel, var::Symbol, class::Symbol;
                       untrans::Bool = false,
                       fourquarter::Bool = false)
    if untrans && fourquarter
        error("Only one of untrans or fourquarter can be true")
    end

    if class in [:obs, :pseudo]
        dict = if class == :obs
            m.observable_mappings
        elseif class == :pseudo
            m.pseudo_observable_mappings
        end
        transform = dict[var].rev_transform

        if transform in [loggrowthtopct_annualized_percapita, loggrowthtopct_annualized, loggrowthtopct, loggrowthtopct_percapita]
            if untrans
                return "Q/Q Log Growth Rate"
            elseif fourquarter
                return "Percent 4Q Growth"
            else
                return "Percent Q/Q Annualized"
            end
        elseif transform in [logleveltopct_annualized_percapita, logleveltopct_annualized]
            if untrans
                return "Log Level"
            elseif fourquarter
                return "Percent 4Q Growth"
            else
                return "Percent Q/Q Annualized"
            end
        elseif transform == quartertoannual
            if untrans
                return "Percent Q/Q"
            else
                return "Percent Annualized"
            end
        elseif transform == identity
            ""
        else
            error("series_ylabel not implemented for transform: $transform")
        end
    elseif class == :stdshocks
        return "Standard Deviations"
    elseif class in [:states, :shocks]
        return ""
    else
        error("Invalid class: " * string(class))
    end
end

function save_plot(p::Plots.Plot, output_file::String = ""; verbose::Symbol = :low)
    if !isempty(output_file)
        output_dir = dirname(output_file)
        !isdir(output_dir) && mkpath(output_dir)
        Plots.savefig(output_file)

        if VERBOSITY[verbose] >= VERBOSITY[:low]
            println("Saved $output_file")
        end
    end
end