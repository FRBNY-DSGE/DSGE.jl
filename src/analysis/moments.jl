"""
```
function load_posterior_moments(m; load_bands = true, include_fixed = false)
```

Load posterior moments (mean, std) of parameters for a particular sample, and optionally
also load 5% and 95% lower and upper bands.

### Keyword Arguments
- `cloud::ParticleCloud`: Optionally pass in a cloud that you want to load the sample from. If the cloud is non-empty then the model object will only be used to find fixed indices and parameter tex labels
- `load_bands::Bool`: Optionally include the 5% and 95% percentiles for the sample of parameters in the returned df
- `include_fixed::Bool`: Optionally include the fixed parameters in the returned df
- `excl_list::Vector{Symbol}`: List parameters by their key that you want to exclude from
loading

### Outputs
- `df`: A dataframe containing the aforementioned moments/bands
"""
function load_posterior_moments(m::AbstractDSGEModel;
                                cloud::Union{SMC.Cloud,DSGE.Cloud,ParticleCloud} = ParticleCloud(m, 0),
                                load_bands::Bool = true,
                                include_fixed::Bool = false,
                                excl_list::Vector{Symbol} = Vector{Symbol}(undef, 0),
                                weighted::Bool = true)
    parameters = m.parameters

    # Read in Posterior Draws
    if cloud_isempty(cloud)
        if get_setting(m, :sampling_method) == :MH
            params = load_draws(m, :full)
            params = get_setting(m, :sampling_method) == :MH ? thin_mh_draws(m, params) : params # TODO
            params = params'
            weights = Vector{Float64}(undef, 0)
            weighted = false
        elseif get_setting(m, :sampling_method) == :SMC
            cloud = get_cloud(m)
            params = get_vals(cloud)
            weights = get_weights(cloud)
        else
            @error "Invalid sampling method"
        end
    else
        params = get_vals(cloud)
        weights = get_weights(cloud)
   end

    # Index out the fixed parameters
    if include_fixed
        free_indices   = 1:n_parameters(m)
    else
        free_indices   = findall(x -> x.fixed == false, parameters)
        parameters = parameters[free_indices]
        params = params[free_indices, :]
    end

    # Remove excluded parameters
    included_indices = setdiff(1:length(free_indices), calculate_excluded_indices(m, excl_list, include_fixed = include_fixed))
    parameters = parameters[included_indices]
    params = params[included_indices, :]
    tex_labels = [DSGE.detexify(parameters[i].tex_label) for i in 1:length(parameters)]

    load_posterior_moments(params, weights, tex_labels, weighted = weighted, load_bands = load_bands)
end

# Aggregating many ParticleClouds to calculate moments across
# multiple SMC estimations
function load_posterior_moments(m::AbstractDSGEModel,
                                clouds::Union{Vector{ParticleCloud},Vector{SMC.Cloud}};
                                load_bands::Bool = true,
                                include_fixed::Bool = false,
                                excl_list::Vector{Symbol} = Vector{Symbol}(undef, 0),
                                weighted::Bool = true)
    n_params = n_parameters(m)
    n_particles = length(clouds[1])
    n_clouds = length(clouds)
    n_free_params = length(Base.filter(x -> x.fixed == false, m.parameters))

    p_mean = Array{Float64}(undef, n_free_params, n_clouds)
    p_lb = Array{Float64}(undef, n_free_params, n_clouds)
    p_ub = Array{Float64}(undef, n_free_params, n_clouds)

    for (i, c) in enumerate(clouds)
        df_i = load_posterior_moments(m; weighted = weighted, cloud = c, load_bands = load_bands, include_fixed = include_fixed, excl_list = excl_list)
        p_mean[:, i] = df_i[!,:post_mean]
        p_lb[:, i] = df_i[!,:post_lb]
        p_ub[:, i] = df_i[!,:post_ub]
    end
    param = load_posterior_moments(m; cloud = clouds[1], load_bands = true, include_fixed = false)[!,:param]
    df = DataFrame(param = param, post_mean = dropdims(mean(p_mean, dims = 2), dims = 2), post_std = dropdims(std(p_mean, dims = 2),dims = 2), post_lb = dropdims(mean(p_lb, dims = 2), dims = 2), post_ub = dropdims(mean(p_ub, dims = 2), dims = 2))

    return df
end

# Base method
# Assumes n_parameters x n_draws
function load_posterior_moments(params::Matrix{Float64}, weights::Vector{Float64}, tex_labels::Vector{String}; weighted::Bool = true, load_bands::Bool = true)
    if size(params, 1) > size(params, 2)
        @warn "`params` argument to load_posterior_moments seems to be oriented incorrectly.
        The argument should be n_params x n_draws. Currently, size(params) = ($(size(params, 1)), $(size(params, 2)))"
    end
    if weighted
        params_mean = vec(mean(params, Weights(weights), dims = 2))
        params_std = vec(std(params, Weights(weights), 2, corrected = false))
    else
        params_mean = vec(mean(params, dims = 2))
        params_std  = vec(std(params, dims = 2))
    end

    df = DataFrame()
    df[!,:param] = tex_labels
    df[!,:post_mean] = params_mean
    df[!,:post_std]  = params_std

    if load_bands
        post_lb = Vector{Float64}(undef, length(params_mean))
        post_ub = similar(post_lb)
        for i in 1:length(params_mean)
            if weighted
                post_lb[i] = quantile(params[i, :], Weights(weights), .05)
                post_ub[i] = quantile(params[i, :], Weights(weights), .95)
            else
                post_lb[i] = quantile(params[i, :], .05)
                post_ub[i] = quantile(params[i, :], .95)
            end
        end

        df[!,:post_lb] = post_lb
        df[!,:post_ub] = post_ub
    end

    return df
end

# For calculating the indices to exclude when making posterior interval plots
function calculate_excluded_indices(m::AbstractDSGEModel, excl_list::Vector{Symbol}; include_fixed::Bool = false)
    parameters = deepcopy(m.parameters)

    # If fixed parameters are excluded, the indices need to be calculated
    # with respect to only the free parameters
    if include_fixed
        return findall(x -> x in excl_list, [parameters[i].key for i in 1:length(parameters)])
    else
        # Remove the fixed parameters
        Base.filter!(x -> x.fixed == false, parameters)
        return findall(x -> x in excl_list, [parameters[i].key for i in 1:length(parameters)])
    end
end

"""
```
moment_tables(m; percent = 0.90, subset_inds = 1:0, subset_string = "",
    groupings = Dict{String, Vector{Parameter}}(), use_mode = false,
    tables = [:prior_posterior_means, :moments, :prior, :posterior],
    caption = true, outdir = "", verbose = :none)
```

Computes prior and posterior parameter moments. Tabulates prior mean, posterior
mean, and bands in various LaTeX tables. These tables will be saved in `outdir`
if it is nonempty, or else in `tablespath(m, \"estimate\")`.

### Inputs

- `m::AbstractDSGEModel`: model object

### Keyword Arguments

- `percent::AbstractFloat`: the percentage of the mass of draws from
  Metropolis-Hastings included between the bands displayed in output tables.
- `subset_inds::AbstractRange{Int64}`: indices specifying the draws we want to use
- `subset_string::String`: short string identifying the subset to be
  appended to the output filenames. If `subset_inds` is nonempty but
  `subset_string` is empty, an error is thrown
- `groupings::Dict{String, Vector{Parameter}}`: see `?parameter_groupings`
- `use_mode::Bool`: use the modal parameters instead of the mean in the
  prior_posterior_means table
- `tables::Vector{Symbol}`: which tables to produce
- `caption::Bool`: whether to include table captions
- `outdir::String`: where to save output tables
- `verbose::Symbol`: desired frequency of function progress messages printed to
  standard out. One of `:none`, `:low`, or `:high`
"""
function moment_tables(m::AbstractDSGEModel; percent::AbstractFloat = 0.90,
                       subset_inds::AbstractRange{Int64} = 1:0, subset_string::String = "",
                       groupings::AbstractDict{String, Vector{Parameter}} = Dict{String, Vector{Parameter}}(),
                       tables = [:prior_posterior_means, :moments, :prior, :posterior],
                       caption = true, outdir = "",
                       verbose::Symbol = :low, use_mode::Bool = false)

    ### 1. Load parameter draws from Metropolis-Hastings

    params = if !isempty(subset_inds)
        # Use subset of draws
        if isempty(subset_string)
            error("Must supply a nonempty subset_string if subset_inds is nonempty")
        end
        load_draws(m, :subset; subset_inds = subset_inds, verbose = verbose)
    else
        # Use all draws
        load_draws(m, :full; verbose = verbose)
    end

    ### 2. Compute posterior moments

    if use_mode
        post_mode = h5read(get_forecast_input_file(m, :mode), "params")
    end
    post_means = vec(mean(params, dims = 1))

    # Save posterior means
    basename = "paramsmean"
    if !isempty(subset_string)
        basename *= "_sub=$(subset_string)"
    end
    filename = workpath(m, "estimate", "$basename.h5")
    h5open(filename, "w") do file
        write(file, "post_means", post_means)
    end
    post_bands = permutedims(find_density_bands(params, percent; minimize = true))

    ### 3. Produce TeX tables

    if :prior_posterior_means in tables
        prior_posterior_table(m, use_mode ? post_mode : post_means;
                              subset_string = subset_string, groupings = groupings,
                              use_mode = use_mode, caption = caption, outdir = outdir)
    end

    if :moments in tables
        prior_posterior_moments_table(m, post_means, post_bands; percent = percent,
                                      subset_string = subset_string, groupings = groupings,
                                      caption = caption, outdir = outdir)
    end

    if :prior in tables
        prior_table(m, groupings = groupings, caption = caption, outdir = outdir)
    end

    if :posterior in tables
        posterior_table(m, post_means, post_bands, percent = percent,
                        subset_string = subset_string, groupings = groupings,
                        caption = caption, outdir = outdir)
    end
    #=
    if :mean_mode_moments in tables
        mean_mode_moments_table(m, post_means, post_bands; percent = percent,
                                subset_string = subset_string, groupings = groupings,
                                caption = caption, outdir = outdir)
    end
    =#
    println(verbose, :low, "Tables are saved as " * tablespath(m, "estimate", "*.tex"))
end

"""
```
moments(θ::Parameter)
```

If θ's prior is a `RootInverseGamma`, τ and ν. Otherwise, returns the mean
and standard deviation of the prior. If θ is fixed, returns `(θ.value, 0.0)`.
"""
function moments(θ::Parameter)
    if θ.fixed
        return θ.value, 0.0
    else
        prior = get(θ.prior)
        if isa(prior, RootInverseGamma)
            return prior.τ, prior.ν
        else
            return mean(prior), std(prior)
        end
    end
end

"""
```
prior_table(m; subset_string = "", groupings = Dict{String, Vector{Parameter}}(),
    caption = true, outdir = "")
```
"""
function prior_table(m::AbstractDSGEModel; subset_string::String = "",
                     caption = true, outdir = "",
             groupings::AbstractDict{String, Vector{Parameter}} = Dict{String, Vector{Parameter}}())

    if isempty(groupings)
        sorted_parameters = sort(m.parameters, by = (x -> x.key))
        groupings[""] = sorted_parameters
    end

    # Open the TeX file
    basename = "priors"
    if !isempty(subset_string)
        basename *= "_sub=$(subset_string)"
    end
    outfile = tablespath(m, "estimate", "$basename.tex")
    if !isempty(outdir)
        outfile = replace(outfile, dirname(outfile) => outdir)
    end
    fid = open(outfile, "w")

    # Write header
    write_table_preamble(fid)

    @printf fid "\\renewcommand*\\footnoterule{}"
    @printf fid "\\vspace*{.5cm}\n"
    @printf fid "{\\small\n"
    @printf fid "\\begin{longtable}{rlrr@{\\hspace{1in}}rlrr}\n"
    if caption
        @printf fid "\\caption{Priors}\n"
    end
    @printf fid "\\label{tab:param-priors}\n"
    @printf fid "\\\\ \\hline\n"

    @printf fid "& Dist & Mean & Std Dev & & Dist & Mean & Std Dev \\\\ \\hline\n"
    @printf fid "\\endhead\n"

    @printf fid "\\hline \\\\\n"
    @printf fid "\\multicolumn{8}{c}{\\footnotesize Note: For Inverse Gamma prior mean and SD, \$\\tau\$ and \$\\nu\$ reported.}\n"
    @printf fid "\\endfoot\n"

    # Map prior distributions to identifying strings
    distid(::Distributions.Uniform) = "Uniform"
    distid(::Distributions.Beta)    = "Beta"
    distid(::Distributions.Gamma)   = "Gamma"
    distid(::Distributions.Normal)  = "Normal"
    distid(::RootInverseGamma)      = "InvG"

    # Write priors
    for group_desc in keys(groupings)
        params = groupings[group_desc]

        # Take out anticipated shock SDs 2 to k - these priors are all the same
        antshock_params = [m[k] for k in [Symbol("σ_r_m$i") for i = 2:n_mon_anticipated_shocks(m)]]
        params = setdiff(params, antshock_params)

        n_params = length(params)
        n_rows = convert(Int, ceil(n_params/2))

        # Write grouping description if not empty
        if !isempty(group_desc)
            @printf fid "\\multicolumn{8}{l}{\\textit{%s}} \\\\[3pt]\n" group_desc
        end

        # Write footnote about standard deviations of anticipated policy shocks
        function anticipated_shock_footnote(θ::Parameter)
            if n_mon_anticipated_shocks(m) > 0 && θ.key == :σ_r_m1
                nantpad          = n_mon_anticipated_shocks_padding(m)
                all_sigmas       = [m[Symbol("σ_r_m$i")]::Parameter for i = 1:nantpad]
                nonzero_sigmas   = Base.filter(θ -> !(θ.fixed && θ.value == 0), all_sigmas)
                n_nonzero_sigmas = length(nonzero_sigmas)

                if n_nonzero_sigmas > 1
                    text = "\$\\sigma_{ant1}\$ through \$\\sigma_{ant$(n_nonzero_sigmas)}\$ all have the same distribution."
                    @printf fid "\\let\\thefootnote\\relax\\footnote{\\centering %s}" text
                end
            end
        end

        for i = 1:n_rows
            # Write left column
            θ = params[i]
            (prior_mean, prior_std) = moments(θ)
            @printf fid "\$%s\$ &" θ.tex_label
            @printf fid " %s &" (θ.fixed ? "-" : distid(get(θ.prior)))
            @printf fid " %0.2f &" prior_mean
            if θ.fixed
                @printf fid " \\scriptsize{fixed} &"
            else
                @printf fid " %0.2f &" prior_std
            end
            anticipated_shock_footnote(θ)

            # Write right column if it exists
            if n_rows + i <= n_params
                θ = params[n_rows + i]
                (prior_mean, prior_std) = moments(θ)
                @printf fid " \$%s\$ &" θ.tex_label
                @printf fid " %s &" (θ.fixed ? "-" : distid(get(θ.prior)))
                @printf fid " %0.2f &" prior_mean
                if θ.fixed
                    @printf fid " \\scriptsize{fixed}"
                else
                    @printf fid " %0.2f" prior_std
                end
                anticipated_shock_footnote(θ)
            else
                @printf fid  "& & &"
            end

            # Add padding after last row in a grouping
            if i == n_rows
                @printf fid " \\\\[3pt]\n"
            else
                @printf fid " \\\\\n"
            end
        end
    end

    # Write footer
    write_table_postamble(fid; small = true)

    # Close file
    close(fid)
end

"""
```
posterior_table(m, post_means, post_bands; percent = 0.9, subset_string = "",
    groupings = Dict{String, Vector{Parameter}}(), caption = true, outdir = "")
```
"""
function posterior_table(m::AbstractDSGEModel, post_means::Vector, post_bands::Matrix;
                         percent::AbstractFloat = 0.9,
                         subset_string::String = "",
                         groupings::AbstractDict{String, Vector{Parameter}} = Dict{String, Vector{Parameter}}(),
                         caption::Bool = true,
                         outdir::String = "")

    if isempty(groupings)
        sorted_parameters = sort(m.parameters, by = (x -> x.key))
        groupings[""] = sorted_parameters
    end

    # Open the TeX file
    basename = "posterior"
    if !isempty(subset_string)
        basename *= "_sub=$(subset_string)"
    end
    outfile = tablespath(m, "estimate", "$basename.tex")
    if !isempty(outdir)
        outfile = replace(outfile, dirname(outfile) => outdir)
    end
    fid = open(outfile, "w")

    # Write header
    write_table_preamble(fid)

    @printf fid "\\renewcommand*\\footnoterule{}"
    @printf fid "\\vspace*{.5cm}\n"
    @printf fid "{\\small\n"
    @printf fid "\\begin{longtable}{rrc@{\\hspace{1in}}rrc}\n"
    if caption
        @printf fid "\\caption{Posteriors}\n"
    end
    @printf fid "\\label{tab:param-posteriors}\n"
    @printf fid "\\\\ \\hline\n"

    lb = (1 - percent)/2 * 100
    ub = 100 - lb
    @printf fid "& Mean & (p%0.0f, p%0.0f) & & Mean & (p%0.0f, p%0.0f) \\\\ \\hline\n" lb ub lb ub
    @printf fid "\\endhead\n"

    @printf fid "\\hline \\\\\n"
    @printf fid "\\endfoot\n"

    # Write priors
    for group_desc in keys(groupings)
        params = groupings[group_desc]
        n_params = length(params)
        n_rows = convert(Int, ceil(n_params/2))

        # Write grouping description if not empty
        if !isempty(group_desc)
            @printf fid "\\multicolumn{6}{l}{\\textit{%s}} \\\\[3pt]\n" group_desc
        end

        for i = 1:n_rows
            # Write left column
            θ = params[i]
            j = m.keys[θ.key]
            @printf fid "\$%s\$ &" θ.tex_label
            @printf fid " %0.2f &" post_means[j]
            if θ.fixed
                @printf fid " \\scriptsize{fixed} &"
            else
                @printf fid " (%0.2f, %0.2f) &" post_bands[j, :]...
            end

            # Write right column if it exists
            if n_rows + i <= n_params
                θ = params[n_rows + i]
                j = m.keys[θ.key]
                (prior_mean, prior_std) = moments(θ)
                @printf fid " \$%s\$ &" θ.tex_label
                @printf fid " %0.2f &" post_means[j]
                if θ.fixed
                    @printf fid " \\scriptsize{fixed}"
                else
                    @printf fid " (%0.2f, %0.2f)" post_bands[j, :]...
                end
            else
                @printf fid " & &"
            end

            # Add padding after last row in a grouping
            if i == n_rows
                @printf fid " \\\\[3pt]\n"
            else
                @printf fid " \\\\\n"
            end
        end
    end

    # Write footer
    write_table_postamble(fid; small = true)

    # Close file
    close(fid)
end

"""
```
prior_posterior_moments_table(m, post_means, post_bands; percent = 0.9,
    subset_string = "", groupings = Dict{String, Vector{Parameter}}(),
    caption = true, outdir = "")
```

Produces a table of prior means, prior standard deviations, posterior means, and
90% bands for posterior draws.
"""
function prior_posterior_moments_table(m::AbstractDSGEModel,
                                       post_means::Vector, post_bands::Matrix;
                                       percent::AbstractFloat = 0.9,
                                       subset_string::String = "",
                                       caption = true, outdir = "",
                 groupings::AbstractDict{String, Vector{Parameter}} = Dict{String, Vector{Parameter}}())

    if isempty(groupings)
        sorted_parameters = sort(m.parameters, by = (x -> x.key))
        groupings[""] = sorted_parameters
    end

    # Open the TeX file
    basename = "moments"
    if !isempty(subset_string)
        basename *= "_sub=$(subset_string)"
    end
    outfile = tablespath(m, "estimate", "$basename.tex")
    if !isempty(outdir)
        outfile = replace(outfile, dirname(outfile) => outdir)
    end
    fid = open(outfile, "w")

    # Write header
    write_table_preamble(fid)

    @printf fid "\\vspace*{.5cm}\n"
    @printf fid "{\\small\n"
    @printf fid "\\begin{longtable}{lcccccc}\n"
    if caption
        @printf fid "\\caption{Parameter Estimates}\n"
    end
    @printf fid "\\\\ \\hline\n"

    # Two-row column names. First row is multicolumn, where entries (i,str) are `i` columns
    # with content `str`.
    colnames0 = [(1,""), (3,"Prior"),(3,"Posterior")]
    colnames = ["Parameter", "Type", "Mean", "SD", "Mean",
                "$(100*percent)\\% {\\tiny Lower Band}",
                "$(100*percent)\\% {\\tiny Upper Band}"]

    @printf fid "\\multicolumn{%d}{c}{%s}" colnames0[1][1] colnames0[1][2]
    for col in colnames0[2:end]
        @printf fid " & \\multicolumn{%d}{c}{%s}" col[1] col[2]
    end
    @printf fid " \\\\\n"
    @printf fid "%s" colnames[1]
    for col in colnames[2:end]
        @printf fid " & %s" col
    end
    @printf fid " \\\\\n"
    @printf fid "\\cmidrule(lr){1-1} \\cmidrule(lr){2-4} \\cmidrule(lr){5-7}\n"
    @printf fid "\\endhead\n"

    @printf fid "\\hline\n"
    @printf fid "\\\\ \\multicolumn{7}{c}{\\footnotesize Note: For Inverse Gamma (IG) prior mean and SD, \$\\tau\$ and \$\\nu\$ reported.}\n"
    @printf fid "\\endfoot\n"

    # Map prior distributions to identifying strings
    distid(::Distributions.Uniform) = "U"
    distid(::Distributions.Beta)    = "B"
    distid(::Distributions.Gamma)   = "G"
    distid(::Distributions.Normal)  = "N"
    distid(::RootInverseGamma) = "IG"

    # Write parameter moments
    # sorted_parameters = sort(m.parameters, by = (x -> x.key))
    for group_desc in keys(groupings)
        params = groupings[group_desc]

        # Write grouping description if not empty
        if !isempty(group_desc)
            @printf fid "\\multicolumn{7}{c}{\\textit{%s}} \\\\[3pt]\n" group_desc
        end

        for param in params
            index = m.keys[param.key]
            (prior_mean, prior_std) = moments(param)

            @printf fid "\$%4.99s\$ & " param.tex_label
            @printf fid "%s & " (param.fixed ? "-" : distid(get(param.prior)))
            @printf fid "%8.3f & " prior_mean
            @printf fid "%8.3f & " prior_std
            @printf fid "%8.3f & " post_means[index]
            @printf fid "%8.3f & %8.3f \\\\\n" post_bands[index, :]...
        end
    end

    # Write footer
    write_table_postamble(fid; small=true)

    # Close file
    close(fid)
end

"""
```
prior_posterior_table(m, post_values; subset_string = "",
    groupings = Dict{String, Vector{Parameter}}(), use_mode = false,
    caption = true, outdir = "")
```

Produce a table of prior means and posterior means or mode.
"""
function prior_posterior_table(m::AbstractDSGEModel, post_values::Vector;
                               subset_string::String = "",
                               groupings::AbstractDict{String,Vector{Parameter}} = Dict{String, Vector{Parameter}}(),
                               caption = true, outdir = "",
                               use_mode::Bool = false)

    if isempty(groupings)
        sorted_parameters = sort(m.parameters, by = (x -> x.key))
        groupings[""] = sorted_parameters
    end

    # Open the TeX file
    basename = use_mode ? "prior_posterior_mode" : "prior_posterior_means"
    if !isempty(subset_string)
        basename *= "_sub=$(subset_string)"
    end
    table_out = tablespath(m, "estimate", "$basename.tex")
    if !isempty(outdir)
        outfile = replace(outfile, dirname(outfile) => outdir)
    end
    fid = open(table_out, "w")

    # Write header
    write_table_preamble(fid)
    @printf fid "\\vspace*{.5cm}\n"
    @printf fid "\\begin{longtable}{ccc}\n"
    if caption
        if use_mode
            @printf fid "\\caption{Parameter Estimates: Prior Mean and Posterior Mode}\n"
        else
            @printf fid "\\caption{Parameter Estimates: Prior and Posterior Means}\n"
        end
    end
    @printf fid "\\\\ \\hline\n"
    @printf fid "Parameter & Prior & Posterior\n"
    @printf fid "\\\\ \\hline\n"
    @printf fid "\\endhead\n"
    @printf fid "\\hline\n"
    @printf fid "\\endfoot\n"

    # Write results
    # sorted_parameters = sort(m.parameters, by = (x -> x.key))

    for group_desc in keys(groupings)
        params = groupings[group_desc]

        # Write grouping description if not empty
        if !isempty(group_desc)
            @printf fid "\\multicolumn{7}{c}{\\textit{%s}} \\\\[3pt]\n" group_desc
        end

        for param in params
            index = m.keys[param.key]

            post_value = if param.fixed
                param.value
            else
                prior = get(param.prior)
                isa(prior, RootInverseGamma) ? prior.τ : mean(prior)
            end

            @printf fid "\$ %4.99s\$ & " param.tex_label
            @printf fid "%8.3f & " post_value
            @printf fid "\\%8.3f \\\\\n" post_values[index]
        end
    end

    # Write footer
    write_table_postamble(fid)

    # Close file
    close(fid)
end

"""
```
find_density_bands(draws::Matrix, percent::AbstractFloat; minimize::Bool=true)
```

Returns a `2` x `cols(draws)` matrix `bands` such that `percent` of the mass of `draws[:,i]`
is above `bands[1,i]` and below `bands[2,i]`.

### Arguments

- `draws`: `ndraws` by `nperiods` matrix of parameter draws (from Metropolis-Hastings, for example)
- `percent`: percent of data within bands (e.g. .9 to get 90% of mass within bands)

### Optional Arguments

- `minimize`: if `true`, choose shortest interval, otherwise just chop off lowest and
  highest (percent/2)
"""
function find_density_bands(draws::AbstractArray, percent::T;
                            minimize::Bool = true) where {T<:AbstractFloat}
    if !(0 <= percent <= 1)
        error("percent must be between 0 and 1")
    end

    ndraws, nperiods = size(draws)

    if ndraws == 1
        band = repeat(draws, outer=(2, 1))
        return band
    end

    band = zeros(2, nperiods)
    n_in_band  = round(Int, percent * ndraws)  # number of draws in the band

    for i in 1:nperiods

        # Sort response for parameter i such that 1st element is largest
        draw_variable_i = draws[:,i]
        sort!(draw_variable_i, rev=true)

        # Search to find the shortest interval containing `percent` of
        # the mass `low` is the index of the largest draw in the band
        # (but the first index to take in `draw_variable_i`, `high` is
        # the smallest (but the highest index to take)

        low = if minimize

            low         = 1
            done        = 0
            j           = 2
            minwidth = draw_variable_i[1] - draw_variable_i[n_in_band]

            while j <= (ndraws - n_in_band + 1)

                newwidth = draw_variable_i[j] - draw_variable_i[j + n_in_band - 1]

                if newwidth < minwidth
                    low = j
                    minwidth = newwidth
                end

                j += 1
            end

            low
        else
            # Chop off lowest and highest percent/2
            ndraws - n_in_band - round(Int, floor(.5*(ndraws-n_in_band)))
        end

        high = low + n_in_band - 1

        if any(ismissing.(draw_variable_i)) || isnan(draw_variable_i[low]) || isnan(draw_variable_i[high])
            band[2,i] = NaN
            band[1,i] = NaN
        else
            band[2,i] = draw_variable_i[low]
            band[1,i] = draw_variable_i[high]
        end
    end
    return band
end

"""
```
find_density_bands(draws::Matrix, percents::Vector{T}; minimize::Bool=true) where T<:AbstractFloat
```

Returns a `2` x `cols(draws)` matrix `bands` such that `percent` of the mass of `draws[:,i]`
is above `bands[1,i]` and below `bands[2,i]`.

### Arguments

- `draws`: Matrix of parameter draws (from Metropolis-Hastings, for example)
- `percent`: percent of data within bands (e.g. .9 to get 90% of mass within bands)

### Optional Arguments

- `minimize`: if `true`, choose shortest interval, otherwise just chop off lowest and
  highest (percent/2)
"""
function find_density_bands(draws::AbstractArray, percents::Vector{T};
                            minimize::Bool = true) where {T<:AbstractFloat}
    bands = DataFrame()

    for p in percents
        out = find_density_bands(draws, p, minimize = minimize)
        bands[!, Symbol("$(100*p)% UB")] = vec(out[2,:])
        bands[!, Symbol("$(100*p)% LB")] = vec(out[1,:])
    end
    bands
end

function write_table_preamble(fid::IOStream)
    @printf fid "\\documentclass[12pt]{article}\n"
    @printf fid "\\usepackage{booktabs}\n"
    @printf fid "\\usepackage[justification=centering]{caption}\n"
    @printf fid "\\usepackage[margin=1in]{geometry}\n"
    @printf fid "\\usepackage{longtable}\n"
    @printf fid "\\usepackage{graphicx}\n"
    @printf fid "\\usepackage{cellspace}\n"
    @printf fid "\\setlength\\cellspacetoplimit{7pt}\n"
    @printf fid "\\setlength\\cellspacebottomlimit{7pt}\n"
    @printf fid "\\begin{document}\n"
    @printf fid "\\pagestyle{empty}\n"
end

# `small`: Whether to print an additional curly bracket after "\end{longtable}" (necessary if
# the table is enclosed by "\small{}")
function write_table_postamble(fid::IOStream; small::Bool=false, tabular::Bool=false)
    if small
        if tabular
            @printf fid "\\end{tabular}}\n"
        else
            @printf fid "\\end{longtable}}\n"
        end
    else
        if tabular
            @printf fid "\\end{tabular}\n"
        else
            @printf fid "\\end{longtable}\n"
        end
    end

    @printf fid "\\end{document}"
end

"""
```
sample_λ(m, pred_dens, θs, T = -1; parallel = true) where S<:AbstractFloat
sample_λ(m, pred_dens, T = -1; parallel = true) where S<:AbstractFloat
```

Computes and samples from the conditional density p(λ_t|θ, I_t, P) for
particle in `θs`, which represents the posterior distribution. The sampled
λ particles represent the posterior distribution p(λ_{t|t} | I_t, P).

If no posterior distribution is passed in, then the function computes
the distribution of λ_{t|t} for a static pool.

### Inputs

- `m::PoolModel{S}`: `PoolModel` object
- `pred_dens::Matrix{S}`: matrix of predictive densities
- `θs::Matrix{S}`: matrix of particles representing posterior distribution of θ
- `T::Int64`: final period for tempered particle filter

where `S<:AbstractFloat`.

### Keyword Argument

- `parallel::Bool`: use parallel computing to compute and sample draws of λ

### Outputs

- `λ_sample::Vector{Float64}`: sample of draws of λs; together with (θ,λ) represents a joint density

```
"""
function sample_λ(m::PoolModel{S}, pred_dens::Matrix{S}, θs::Matrix{S}, T::Int64 = -1;
                  parallel::Bool = false,
                  tuning0::Dict{Symbol,Any} = Dict{Symbol,Any}()) where S<:AbstractFloat
    # Check size and orientation of θs is correct: assume particle_num x parameter_num
    if length(m.parameters) != size(θs,2)
        error("number of parameters in PoolModel do not match number of parameters in matrix of posterior draws of θ")
    end

    # Initialize necessary objects
    if T <= 0
         error("T must be positive") # No period provided or is invalid
    end
    tuning = isempty(tuning0) ? deepcopy(get_setting(m, :tuning)) : deepcopy(tuning0)
    tuning[:get_t_particle_dist] = true
    tuning[:parallel] = false
    tuning[:allout] = false

    # Sample from p(λ|θ, I_t^P, P) for each θ in posterior
    data = (T == 1) ? reshape(pred_dens[:,1], 2, 1) : pred_dens[:,1:T]
    if parallel
        # Send variables to workers to avoid issues with serialization
        # across workers with different Julia system images
        sendto(workers(), m = m)
        sendto(workers(), data = data)
        sendto(workers(), θs = θs)
        sendto(workers(), tuning = tuning)

        # Sample from λ distribution given θ
        λ_sample = @sync @distributed (vcat) for i in 1:size(θs,1)
            sample_λ(m, data, vec(θs[i,:]), tuning)
        end
        if sum(λ_sample) == 0
            error("Sums to zero")
        end
        λ_sample = Array(λ_sample)
    else
        # Same as above but everything is serialized
        λ_sample = Vector{Float64}(undef, size(θs,1))
        for i in 1:size(θs,1)
            λ_sample[i] = sample_λ(m, data, vec(θs[i,:]), tuning)
        end
    end

    return λ_sample
end

# This function actually does the sampling for a given θ,
# but we provide a wrapper for an easier user experience
function sample_λ(m::PoolModel{S}, data::Matrix{S}, θ::Vector{S},
                  tuning::Dict{Symbol,Any}) where S<:AbstractFloat
    update!(m, θ)
    loglik, λ_particles, λ_weights = DSGE.filter(m, data; tuning = tuning)
    return λ_particles[size(data,2)][1,DSGE.sample(DSGE.Weights(λ_weights[:,end]))]
end

function sample_λ(m::PoolModel{S}, pred_dens::Matrix{S}, T::Int64 = -1;
                  parallel::Bool = false,
                  filestring_addl::Vector{String} = Vector{String}(undef,0),
                  tuning0::Dict{Symbol,Any} = Dict{Symbol,Any}()) where S<:AbstractFloat

    # Check m is static
    if get_setting(m, :weight_type) != :static
        error("PoolModel is not static. Set the keyword argument weight_type = :static to create a static PoolModel object.")
    end

    # Initialize necessary objects
    if T <= 0
         error("T must be positive") # No period provided or is invalid
    end
    tuning = isempty(tuning0) ? deepcopy(get_setting(m, :tuning)) : deepcopy(tuning0)
    tuning[:get_t_particle_dist] = true
    tuning[:parallel] = parallel
    tuning[:allout] = false
    orig_samp_method = get_setting(m, :sampling_method)
    m <= Setting(:sampling_method, :MH)

    # Compute posterior from a static pool
    data = (T == 1) ? reshape(pred_dens[:,1], 2, 1) : pred_dens[:,1:T] # make sure it is matrix
    estimate(m, data; filestring_addl = filestring_addl, proposal_covariance = ones(1,1))
    m <= Setting(:sampling_method, orig_samp_method)

    return h5read(rawpath(m, "estimate", "mhsave.h5", filestring_addl), "mhparams")
end

"""
```
compute_Eλ(m, h, λvec, θmat = [], weights = [];
    current_period = true, parallel = true) where T<:AbstractFloat
```

Computes and samples from the conditional density p(λ_t|θ, I_t, P) for
particle in `θs`, which represents the posterior distribution.

### Inputs

- `m::PoolModel{T}`: `PoolModel` object
- `h::Int64`: forecast horizon
- `λvec::Vector{T}`: vector of particles of λ samples from (θ,λ) joint distribution
- `θmat::Matrix{T}': matrix of posterior parameter samples
- `weights::Vector{T}`: weights of λ particles, defaults to equal weights

### Keyword Argument

- `current_period::Bool`: compute Eλ for current period t
- `parallel::Bool`: use parallel computing to compute and sample λ
- `get_dpp_pred_dens::Bool`: compute predictive densities according to dynamic prediction pools

### Outputs

- `λhat_tplush::Float64`: E[λ_{t+h|t} | I_t^P, P]
- `λhat_t::Float64`: E[λ_{t|t} | I_t^P, P]
```
"""
function compute_Eλ(m::PoolModel{T}, h::Int64, λvec::Vector{T},
                    θmat::Matrix{T} = Matrix{Float64}(undef,0,0),
                    weights::Vector{T} = Vector{Float64}(undef,0);
                    current_period::Bool = true, parallel::Bool = false) where T<:AbstractFloat

    # Set up
    if isempty(weights)
        weights = ones(length(λvec)) # assume equal weights
    end
    update_param = !isempty(θmat)
    if update_param # save so we can reset the parameters of the PoolModel
        old_param = map(x -> x.value, m.parameters)
    end
    λhat_t = if current_period mean(λvec .* weights) end # compute expected lambda in current period t

    # Push forward states and compute mean
    if parallel
        # Send data to all workers
        sendto(workers(), m = m)
        sendto(workers(), h = h)
        if update_param
            sendto(workers(), θmat = θmat)

            # Propagate forward!
            λ_vec = @sync @distributed (vcat) for i in 1:length(λvec)
                propagate_λ(λvec[i], h, m, vec(θmat[i,:]))
            end
        else
            λ_vec = @sync @distributed (vcat) for i in 1:length(λvec)
                propagate_λ(λvec[i], h, m)
            end
        end
    else
        if update_param
            λ_vec = map(i -> propagate_λ(λvec[i], h, m, vec(θmat[i,:])), collect(1:length(λvec)))
        else
            λ_vec = map(x -> propagate_λ(x, h, m), λvec)
        end
    end
    λhat_tplush = mean(λ_vec .* weights)

    if update_param
        # Update m to take original parameters
        DSGE.update!(m, old_param)
    end

    if current_period
        return λhat_tplush, λhat_t
    else
        return λhat_tplush
    end
end

"""
```
propgate_λ(λvec, h, m, θvec) where T<:AbstractFloat
```

Propagates a λ particle h periods forward.

### Inputs

- `λ::T`: λ sample from (θ,λ) joint distribution
- `h::Int64`: forecast horizon
- `m::PoolModel`: PoolModel object
- `θvec::Vector{T}`: optional vector of parameters to update PoolModel

```
"""
function propagate_λ(λ::T, h::Int64, m::PoolModel,
                     θvec = Vector{Float64}(undef,0)) where T<:AbstractFloat
    if !isempty(θvec)
        update!(m, θvec)
    end
    Φ, ~, ~ = solve(m)
    for j in 1:h
        λ = Φ([λ; 1 - λ], [0.])[1]
    end
    return λ
end
