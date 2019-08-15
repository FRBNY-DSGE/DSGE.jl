"""
```
packet_help()
```

Print sample packet-writing code.
"""
function packet_help()
    standard_code =
        """
          fid = open(\"path/to/tex/file.tex\", \"w\")
          write_preamble(fid, \"Packet Title\", \"Packet Author\",
                         addl_packages = [\"package1\", \"package2\"])
          write_spec_section(fid, m, purpose = \"Packet Purpose\")
          write_estimation_section(fid, m)
          write_forecast_section(fid, m, input_type, cond_type)
          write_postamble(fid)
          close(fid)
        """
    println("Sample standard packet code:")
    println()
    println(standard_code)
    println()

    comparison_code =
        """
          fid = open(joinpath(\"path/to/tex/file.tex\"), \"w\")
          write_preamble(fid, \"Packet Title\", \"Packet Author\",
                         addl_packages = [\"package1\", \"package2\"])
          write_spec_section(fid, m1, m2, purpose = \"Packet Purpose\")
          write_forecast_section(fid, m1, input_type1, cond_type1, m2, input_type2, cond_type2,
                                 forecast_plotroot = \"path/to/forecast/comparison/plots\",
                                 shockdec_plotroot1 = \"path/to/model1/shockdec/plots\",
                                 shockdec_plotroot2 = \"path/to/model2/shockdec/plots\")
          write_postamble(fid)
          close(fid)
        """
    println("Sample comparison packet code:")
    println()
    println(comparison_code)
end

"""
```
write_preamble(fid, title = "", author = ""; addl_packages = [])
```

Begin LaTeX document, load packages, and make title.
"""
function write_preamble(fid::IOStream, title::String = "", author::String = "";
                        addl_packages::Vector{String} = String[])
    @printf fid "\\documentclass{article}\n"
    @printf fid "\n"

    @printf fid "\\usepackage{booktabs} %% For \\midrule in parameter tables\n"
    @printf fid "\\usepackage{fullpage} %% Set all margins to 1 inch\n"
    @printf fid "\\usepackage{graphicx}\n"
    @printf fid "\\usepackage[bookmarksopen=true]{hyperref} %% Create section bookmarks\n"
    @printf fid "\\usepackage{longtable} %% Break tables across pages\n"
    @printf fid "\\usepackage{standalone} %% Skips extra preambles in inputted files\n"
    @printf fid "\\usepackage{tabularx} %% Force table to page width\n"
    @printf fid "\\usepackage{url} %% For \\path command\n"
    @printf fid "\n"

    for package in addl_packages
        @printf fid "\\usepackage{%s}\n" package
    end
    if !isempty(addl_packages)
        @printf fid "\n"
    end

    @printf fid "\\title{%s}\n" title
    @printf fid "\\author{%s}\n" author
    @printf fid "\n"

    @printf fid "\\begin{document}\n"
    @printf fid "\\maketitle\n"
end

function write_overview_section(fid::IOStream, overview::String)
    @printf fid "\n\n"
    @printf fid "\\section{Overview}\n"
    @printf fid "\n"
    @printf fid "%s\n" overview
end


"""
```
write_postamble(fid)
```

End LaTeX document.
"""
function write_postamble(fid::IOStream)
    @printf fid "\n\n"
    @printf fid "\\end{document}\n"
end

"""
```
month_label(m)
```

Create the correct month labels for a forecast.
"""
function month_label(m::AbstractModel)
    date =
    # If within the same year
    if (Dates.year(date_forecast_start(m)) - 2000) == parse(data_vintage(m)[1:2])
        max(Dates.month(date_forecast_start(m)), parse(data_vintage(m)[3:4]))
    elseif (Dates.year(date_forecast_start(m)) - 2000) < parse(data_vintage(m)[1:2])
        parse(data_vintage(m)[3:4])
    else
        Dates.month(date_forecast_start(m))
    end

    months = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul",
              "Aug", "Sep", "Oct", "Nov", "Dec"]
    monthnums = collect(1:12)
    ind = findin(monthnums, date)[1]
    return months[ind]
end
