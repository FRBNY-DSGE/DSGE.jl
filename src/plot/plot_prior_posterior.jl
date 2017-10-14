"""
```
plot_prior_posterior(m, param_key; kwargs...)

plot_prior_posterior(m, param_keys = map(x -> x.key, m.parameters),
    include_fixed = false, verbose = :low, kwargs...)
```

Plot prior distribution and histogram of posterior draws for specified
`param_key` or `param_keys`.

### Inputs

- `m::AbstractModel`
- `param_key::Symbol` or `param_keys::Vector{Symbol}`: e.g. `:α` or `[:α, :ζ_p]`

### Keyword Arguments

- `verbose::Symbol`

See `?priorpost` for additional keyword arguments, all of which can be passed
into `plot_prior_posterior`.

**Method 2 only:**

- `include_fixed::Bool`: whether to plot fixed parameters
"""
function plot_prior_posterior(m::AbstractModel, param_key::Symbol; kwargs...)

    plots = plot_prior_posterior(m, [param_key]; include_fixed = true, kwargs...)
    return plots[param_key]
end

function plot_prior_posterior(m::AbstractModel,
                              param_keys::Vector{Symbol} = map(x -> x.key, m.parameters);
                              include_fixed::Bool = false,
                              verbose::Symbol = :low,
                              kwargs...)
    # Load parameter draws from Metropolis-Hastings
    posterior_draws = load_draws(m, :full)

    # Loop through parameters
    plots = OrderedDict{Symbol, Plots.Plot}()
    for key in param_keys
        if !include_fixed && param.fixed
            continue
        end

        i = findfirst(x -> x.key == key, m.parameters)
        if i > 0
            param = m.parameters[i]
            posterior = posterior_draws[:, i]

            # Call recipe
            plots[key] = priorpost(param, posterior; verbose = verbose, kwargs...)

            # Save plot
            output_file = figurespath(m, "estimate", "prior_posterior_" * detexify(string(param.key)) * ".pdf")
            save_plot(plots[key], output_file, verbose = verbose)

        else
            error("Parameter not found: " * string(param_key))
        end
    end
    return plots
end