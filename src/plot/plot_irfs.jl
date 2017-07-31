"""
```
plot_irfs(shock, var, mb; bands = [\"90.0\"], output_file = "", title = "")

plot_irfs(shock, vars, mb; bands = [\"90.0%\"], output_file = "", titles = [],
    layout = Plots.EmptyLayout())
```

Plots the responses of either a single `var` (method 1) or all provided `vars`
(method 2) to `shock`, given a `MeansBands` object `mb`. By default, only 90%
bands are plotted. In method 2, a single plot with subplots for each variable is
returned and/or saved.

### Inputs

- `shock::Symbol`: e.g. `:g_sh`
- `var::Symbol` or `vars::Vector{Symbol}`: response variable(s), e.g. `:obs_gdp`
- `mb::MeansBands`

### Keyword Arguments

- `bands::Vector{String}`
- `output_file::String`: if specified, plot will be saved there as a PDF
- `title::String` or `titles::Vector{String}`
- `layout::Plots.AbstractLayout` (method 2 only)

### Output

- `p::Plot`
"""
function plot_irfs(shock::Symbol, var::Symbol, mb::MeansBands;
                   bands::Vector{String} = ["90.0%"],
                   output_file::String = "",
                   title::String = "")

    # Plot title
    if isempty(title)
        title = string(var)
    end
    p = plot(title = title, margin = 10px)

    varshock = Symbol("$(var)__$(shock)")

    # Plot 90% bands
    for pct in bands
        plot!(p, mb.bands[varshock][Symbol("$pct UB")],
              fillto = mb.bands[varshock][Symbol("$pct LB")],
              label = "", color = :blue, α = 0.10)
    end

    # Plot mean
    plot!(p, mb.means[varshock], label = "", linewidth = 2, linecolor = :black)

    # Save if `output_file` provided
    save_plot(p, output_file)

    return p
end

function plot_irfs(shock::Symbol, vars::Vector{Symbol}, mb::MeansBands;
                   bands::Vector{String} = ["90.0%"],
                   output_file::String = "",
                   titles::Vector{String} = Vector{String}(),
                   layout::Plots.AbstractLayout = Plots.EmptyLayout())

    # Determine title
    if isempty(titles)
        titles = map(string, vars)
    end

    plots = OrderedDict{Symbol, Plots.Plot}()

    # Iterate through variables
    for (var, title) in zip(vars, titles)
        plots[var] = plot_irfs(shock, var, mb, bands = bands, title = title)
    end

    # Plot all variables together
    p = if isa(layout, Plots.EmptyLayout)
        plot(collect(values(plots))...)
    else
        plot(collect(values(plots))..., layout = layout)
    end

    # Save if `output_file` provided
    save_plot(p, output_file)

    return p
end