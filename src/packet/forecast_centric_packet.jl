"""
```
write_forecast_centric_packet(m, input_type, cond_type,
    output_vars = [:forecastobs, :forecastpseudo, :shockdecobs, :shockdecpseudo];
    sections = [:estimation, :forecast]
    forecast_string = "", outdir = joinpath(saveroot(m), \"Packets\", spec(m), subspec(m)))
```

Write standard estimation and forecast result packet to `outdir`, with an emphasis
on the forecasting results. The main difference with the standard packet is that the
prior posterior plots at placed at the very end.
"""
function write_forecast_centric_packet(m::AbstractModel, input_type::Symbol, cond_type::Symbol,
                               output_vars::Vector{Symbol} = [:forecastobs, :forecastpseudo,
                                                              :shockdecobs, :shockdecpseudo];
                               sections::Vector{Symbol} = [:estimation, :forecast],
                               forecast_string::String = "",
                               outdir::String = joinpath(saveroot(m), "Packets", spec(m), subspec(m)),
                               purpose::String = "")
    @assert issubset(sections, [:estimation, :forecast, :irf]) "Section specified in `section` kwarg is not supported. Must be a subset of [:estimation, :forecast, :irf]."

    # Title and authors
    title = "Estimation and Forecasting Results \\\\ " * DSGE.description(m)
    authors = "User"

    # Compute file name
    cdvt_str = cond_type == :none ? "" : "_cdvt=" * cond_vintage(m)
    fn = joinpath(outdir, "results_cond=" * string(cond_type) * cdvt_str * "_fcid=" * forecast_string * "_vint=" * data_vintage(m) * ".tex")
    isdir(dirname(fn)) || mkpath(dirname(fn))

    # Write packet
    open(fn, "w") do fid
        write_preamble(fid, title, authors)
        write_spec_section(fid, m, purpose = purpose)
        if :estimation in sections
            write_estimation_table_section(fid, m)
        end
        if :forecast in sections
            write_forecast_section(fid, m, input_type, cond_type,
                                   setdiff(output_vars, [:irfstates, :irfobs, :irfpseudo]),
                                   forecast_string = forecast_string)
        end
        if :irf in sections
            write_irf_section(fid, m , input_type, cond_type, output_vars,
                              forecast_string = forecast_string)
        end
        if :estimation in sections
            write_estimation_plot_section(fid, m)
        end
        write_postamble(fid)
    end
    println("Wrote " * fn)

end

"""
```
write_estimation_table_section(fid, m; plotroot = "")
```

Write parameter moment tables.
"""
function write_estimation_table_section(fid::IOStream, m::AbstractModel;
                                  plotroot::String = "")
    @printf fid "\n\n"
    @printf fid "\\clearpage\n"
    @printf fid "\\section{Estimation Moments}\n"
    @printf fid "\n"

    @printf fid "\\input{%s}\n" tablespath(m, "estimate", "moments.tex")


end

"""
```
write_estimation_plot_section(fid, m; plotroot = "")
```

Write prior/posterior plots.
"""
function write_estimation_plot_section(fid::IOStream, m::AbstractModel;
                                  plotroot::String = "")

    @printf fid "\\clearpage\n"
    @printf fid "\\section{Estimation Histograms}\n"
    write_estimation_plots(fid, m, plotroot = plotroot)
    @printf fid "\n"

end
