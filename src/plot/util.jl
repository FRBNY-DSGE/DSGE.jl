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

function describe_series(m::AbstractDSGEModel, var::Symbol, class::Symbol;
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
        replace(string(var), r"_sh$" => "")
    else
        error("Invalid class: " * string(class))
    end

    detexify ? DSGE.detexify(res) : res
end

function describe_series(m::AbstractDSGEVARModel, var::Symbol, class::Symbol;
                         detexify::Bool = false)
    return describe_series(get_dsge(m), var, class; detexify = detexify)
end

function series_ylabel(m::AbstractDSGEModel, var::Symbol, class::Symbol;
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
        elseif transform in [quartertoannual]
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

function save_plot(output_file::String = ""; verbose::Symbol = :low)
    if !isempty(output_file)
        output_dir = dirname(output_file)
        !isdir(output_dir) && mkpath(output_dir)
        Plots.savefig(output_file)

        println(verbose, :low, "Saved $(abspath(output_file))")
    end
end

function save_plot(p::Plots.Plot, output_file::String = ""; verbose::Symbol = :low)
    if !isempty(output_file)
        output_dir = dirname(output_file)
        !isdir(output_dir) && mkpath(output_dir)
        Plots.savefig(p, output_file)

        println(verbose, :low, "Saved $(abspath(output_file))")
    end
end

function rescale_annotations_fontsize!(p::Plots.Plot{Plots.GRBackend}, fontsize::Int64)
    for i in 1:length(p.series_list)
        if typeof(p.series_list[i].plotattributes[:series_annotations]) == Nothing
            continue
        else
            p.series_list[i].plotattributes[:series_annotations].font.pointsize = fontsize
        end
    end
end

# horiz_position can be :hcenter, :left, :right
# vert_position can be :vcenter, :top, :bottom
# reverse_position Bool to indicate if you want left and right to actually correspond to left and right..
function relocate_annotations!(p::Plots.Plot{Plots.GRBackend},
                               horiz_position::Symbol, vert_position::Symbol,
                               reverse_position::Bool = true)
    if reverse_position
        horiz_position =
        if horiz_position == :right
            :left
        elseif horiz_position == :left
            :right
        else
            :hcenter
        end
        vert_position =
        if vert_position == :top
            :bottom
        elseif vert_position == :bottom
            :top
        else
            :vcenter
        end
    end

    for i in 1:length(p.series_list)
        if typeof(p.series_list[i].plotattributes[:series_annotations]) == Nothing
            continue
        else
            p.series_list[i].plotattributes[:series_annotations].font.halign = horiz_position
            p.series_list[i].plotattributes[:series_annotations].font.valign = vert_position
        end
    end
end

function date_to_float(datevec::Vector{Date}, reps::Int64 = 1)
    date_floats = zeros(reps, length(datevec))
    tofloats = Dict{String,Float64}("03" => .25, "06" => .5, "09" => .75, "12" => 1.)
    for j = 1:length(datevec)
        tmp = string(datevec[j])
        s = tmp[6:7]
        num = parse(Float64,tmp[1:4]) + tofloats[s]
        date_floats[:,j] .= num
    end
    return date_floats
end
