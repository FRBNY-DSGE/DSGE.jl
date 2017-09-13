"""
```
plot_prior_posterior(m::AbstractModel; include_fixed = false, sub_dir = "")

plot_prior_posterior(m::AbstractModel, param_key::Symbol, posterior_draws::Matrix{T};
    output_file = "")

plot_prior_posterior(param::Parameter, posterior_draws::Vector{T};
    output_file = "")
```

Plot prior distribution and histogram of posterior draws for either all
parameters (method 1) or the specified `param_key` or `param` (methods 2 and
3). In method 1, `include_fixed` determines whether fixed parameters are also
plotted.
"""
function plot_prior_posterior(m::AbstractModel;
                              include_fixed::Bool = false,
                              sub_dir::String = "")
    # Load parameter draws from Metropolis-Hastings
    posterior_draws = load_draws(m, :full)

    for i in 1:n_parameters(m)
        param = m.parameters[i]

        if !include_fixed && param.fixed
            continue
        end

        posterior = posterior_draws[:, i]
        output_file = isempty(sub_dir) ? figurespath(m, "estimate", "prior_posterior_" * detexify(string(param.key)) * ".pdf") : figurespath(m, "estimate", sub_dir*"/prior_posterior_" * detexify(string(param.key)) * ".pdf")
        plot_prior_posterior(param, posterior, output_file = output_file)
    end
end

function plot_prior_posterior{T<:AbstractFloat}(m::AbstractModel, param_key::Symbol,
                                                posterior_draws::Matrix{T};
                                                output_file::String = "")
    i = findfirst(x -> x.key == param_key, m.parameters)
    if i > 0
        param = m.parameters[i]
        posterior = posterior_draws[:, i]
        output_file = figurespath(m, "estimate", "prior_posterior_" * detexify(string(param.key)) * ".pdf")
        plot_prior_posterior(param, posterior, output_file = output_file)
    else
        error("Parameter not found: " * string(param_key))
    end
end

function plot_prior_posterior{T<:AbstractFloat}(param::Parameter, posterior_draws::Vector{T};
                                                output_file::String = "")
    # Plot posterior
    label = if param.fixed
        "Posterior: " * describe_prior(param)
    else
        "Posterior"
    end
    p = histogram(posterior_draws, normalize = true, label = label, legend = :bottomright)

    # Plot prior
    prior_label = "Prior: " * describe_prior(param)
    if param.fixed
        plot!(p, [Plots.xlims(p)...], [1, 1], label = prior_label, linewidth = 4)
    elseif !param.fixed && !isnull(param.prior)
        prior = get(param.prior)
        plot!(p, typeof(prior), prior, label = prior_label, linewidth = 4)
    end

    # Add title
    title!(p, param.tex_label)

    # Save if output_file provided
    save_plot(p, output_file)

    return p
end
