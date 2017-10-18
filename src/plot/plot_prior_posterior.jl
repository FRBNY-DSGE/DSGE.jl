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
        i = findfirst(x -> x.key == key, m.parameters)
        if i > 0
            param = m.parameters[i]
            posterior = posterior_draws[:, i]

            # Skip fixed parameter if desired
            if !include_fixed && param.fixed
                continue
            end

            # Call recipe
            plots[key] = priorpost(param, posterior; kwargs...)

            # Save plot
            output_file = figurespath(m, "estimate", "prior_posterior_" * detexify(string(param.key)) * ".pdf")
            save_plot(plots[key], output_file, verbose = verbose)

        else
            error("Parameter not found: " * string(param_key))
        end
    end
    return plots
end

@userplot PriorPost

"""
```
priorpost(param, posterior_draws)
```

User recipe called by `plot_prior_posterior`.

### Inputs

- `param::Parameter`
- `posterior_draws::AbstractVector`

### Keyword Arguments

- `prior_color`
- `posterior_color`

Additionally, all Plots attributes (see docs.juliaplots.org/latest/attributes)
are supported as keyword arguments.
"""
priorpost

@recipe function f(pp::PriorPost;
                   prior_color = :red,
                   posterior_color = :blue)
    # Error checking
    if length(pp.args) != 2 || !(typeof(pp.args[1]) <: Parameter) || !(typeof(pp.args[2]) <: AbstractVector)
        error("priorpost must be given a Parameter and an AbstractVector. Got $(typeof(pp.args))")
    end

    param, posterior_draws = pp.args

    title  --> param.tex_label
    legend --> :bottomright

    # Posterior
    @series begin
        seriestype := :histogram
        normalize  := true
        label      := param.fixed ? "Posterior: " * DSGE.describe_prior(param) : "Posterior"
        color      := posterior_color

        posterior_draws
    end

    # Prior
    @series begin
        label      := "Prior: " * DSGE.describe_prior(param)
        color      := prior_color
        linewidth --> 4

        if param.fixed
            @show keys(plotattributes)
            [minimum(posterior_draws), maximum(posterior_draws)], [1, 1]
        elseif !param.fixed && !isnull(param.prior)
            prior = get(param.prior)
            typeof(prior), prior
        else
            error("Parameter must either be fixed or have a non-null prior")
        end
    end
end
