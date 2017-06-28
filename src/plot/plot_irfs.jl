"""
```
plot_irfs(shock, vars, mb; bands = [\"90.0%\"], output_file = "")
```

Plots the responses of `vars` to `shock`, given a `MeansBands` object `mb`. By
default, only 90% bands are plotted.

### Inputs

- `shock::Symbol`: e.g. `:g_sh`
- `vars::Vector{Symbol}`: response variables, e.g. `[:obs_gdp, :obs_corepce]`
- `mb::MeansBands`

### Keyword Arguments

- `bands::Vector{String}`
- `output_file::String`: if `output_file` is provided in the form \"base.ext\",
  the individual plots will also each be saved as \"base__var.ext\"

### Output

- `allplots::OrderedDict{Symbol, Plot}`: indexed by variable name
"""
function plot_irfs(shock::Symbol, vars::Vector{Symbol}, mb::MeansBands;
                   bands::Vector{String} = ["90.0%"],
                   output_file::String = "")
    # Split `output_file` into base and extension
    if !isempty(output_file)
        strs = split(output_file, ".")
        base = join(strs[1:end-1], ".")
        ext  = "." * strs[end]
        save_plots = true
    end

    allplots = OrderedDict{Symbol, Plot}()

    # Iterate through variables
    for var in vars
        p = Plots.plot(title = string(var), margin = 10px)
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
        if save_plots
            this_file = base * "__$var" * ext
            Plots.savefig(this_file)
            println("Saved $this_file")
        end

        allplots[var] = p
    end

    return allplots
end