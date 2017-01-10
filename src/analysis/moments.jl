"""
```
moment_tables(m; percent = 0.90, subset_inds = [], subset_string = "",
    verbose = :none)
```

Computes prior and posterior parameter moments. Tabulates prior mean, posterior mean, and
bands in various LaTeX tables stored `tablespath(m)`.

### Inputs

- `m::AbstractModel`: model object

### Keyword Arguments

- `percent::AbstractFloat`: the percentage of the mass of draws from
  Metropolis-Hastings included between the bands displayed in output tables.
- `subset_inds::Vector{Int}`: indices specifying the draws we want to use
- `subset_string::AbstractString`: short string identifying the subset to be
  appended to the output filenames. If `subset_inds` is nonempty but
  `subset_string` is empty, an error is thrown
- `verbose::Symbol`: desired frequency of function progress messages printed to
  standard out. One of `:none`, `:low`, or `:high`
"""
function moment_tables(m::AbstractModel; percent::AbstractFloat = 0.90,
                       subset_inds::Vector{Int} = Vector{Int}(), subset_string::AbstractString = "",
                       verbose::Symbol = :low)

    ### 1. Collect prior information

    prior_means, prior_std, prior_dist = collect_prior_info(m)


    ### 2. Load parameter draws from Metropolis-Hastings

    filename = get_input_file(m, :full)
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


    ### 3. Compute posterior moments

    post_means = vec(mean(params, 1))
    post_bands = find_density_bands(params, percent; minimize = true)'

    # Save posterior means
    basename = "paramsmean"
    if !isempty(subset_string)
        basename *= "_sub=$(subset_string)"
    end
    filename = workpath(m, "estimate", "$basename.h5")
    h5open(filename, "w") do file
        write(file, "post_means", post_means)
    end


    ### 4. Produce TeX tables

    prior_posterior_moments_table(m, prior_dist, prior_means, prior_std, post_means, post_bands;
                                  percent = percent, subset_string = subset_string)
    prior_posterior_means_table(m, prior_means, post_means; subset_string = subset_string)

    if VERBOSITY[verbose] >= VERBOSITY[:low]
        @printf "Tables are saved as %s.\n" tablespath(m, "estimate", "*.tex")
    end
end

"""
```
collect_prior_info(m)
```

### Inputs

- `m::AbstractModel`: model object

### Outputs

- `prior_means::Vector{AbstractFloat}`
- `prior_std::Vector{AbstractFloat}`
- `prior_dist::Vector{AbstractString}`
"""
function collect_prior_info(m::AbstractModel)
    n_params    = length(m.parameters)
    prior_means = zeros(n_params)
    prior_std   = zeros(n_params)
    prior_dist  = Vector{AbstractString}(n_params)

    # Map prior distributions to identifying strings
    distid(::Distributions.Beta)    = "B"
    distid(::Distributions.Gamma)   = "G"
    distid(::Distributions.Normal)  = "N"
    distid(::DSGE.RootInverseGamma) = "IG"

    for i = 1:n_params
        θ = m.parameters[i]
        if θ.fixed
            # If param is fixed, we simply set mean to its value and std to 0
            prior_means[i] = θ.value
            prior_std[i]   = 0
            prior_dist[i]  = "-"
        else
            prior = get(θ.prior)
            α, β  = moments(prior)
            prior_means[i] = α
            prior_std[i]   = β
            prior_dist[i]  = distid(prior)
        end
    end

    return prior_means, prior_std, prior_dist
end

"""
```
prior_posterior_moments_table(m, prior_dist, prior_means, prior_std, post_means,
    post_bands; percent = 0.9, subset_string = "")
```

Produces a table of prior means, prior standard deviations, posterior means, and
90% bands for posterior draws. Saves to:

```
tablespath(m, \"estimate\", \"moments[_sub=\$subset_string].tex\")
```
"""
function prior_posterior_moments_table(m::AbstractModel,
                                       prior_dist::Vector, prior_means::Vector, prior_std::Vector,
                                       post_means::Vector, post_bands::Matrix;
                                       percent::AbstractFloat = 0.9,
                                       subset_string::AbstractString = "")
    # Collect inputs into one matrix
    outmat = hcat(prior_dist, prior_means, prior_std, post_means, post_bands)

    # Open the TeX file
    basename = "moments"
    if !isempty(subset_string)
        basename *= "_sub=$(subset_string)"
    end
    moments_out = tablespath(m, "estimate", "$basename.tex")
    moments_fid = open(moments_out, "w")

    # Write header
    write_table_preamble(moments_fid)

    @printf moments_fid "\\vspace*{.5cm}\n"
    @printf moments_fid "{\\small\n"
    @printf moments_fid "\\begin{longtable}{lcccccc}\n"
    @printf moments_fid "\\caption{Parameter Estimates}\n"
    @printf moments_fid "\\\\ \\hline\n"

    # Two-row column names. First row is multicolumn, where entries (i,str) are `i` columns
    # with content `str`.
    colnames0 = [(1,""), (3,"Prior"),(3,"Posterior")]
    colnames = ["Parameter", "Type", "Mean", "SD", "Mean",
                "$(100*percent)\\% {\\tiny Lower Band}",
                "$(100*percent)\\% {\\tiny Upper Band}"]

    @printf moments_fid "\\multicolumn{%d}{c}{%s}" colnames0[1][1] colnames0[1][2]
    for col in colnames0[2:end]
        @printf moments_fid " & \\multicolumn{%d}{c}{%s}" col[1] col[2]
    end
    @printf moments_fid " \\\\\n"
    @printf moments_fid "%s" colnames[1]
    for col in colnames[2:end]
        @printf moments_fid " & %s" col
    end
    @printf moments_fid " \\\\\n"
    @printf moments_fid "\\cmidrule(lr){1-1} \\cmidrule(lr){2-4} \\cmidrule(lr){5-7}\n"
    @printf moments_fid "\\endhead\n"

    @printf moments_fid "\\hline\n"
    @printf moments_fid "\\\\ \\multicolumn{7}{c}{\\footnotesize Note: For Inverse Gamma (IG) prior mean and SD, \$\\tau\$ and \$\\nu\$ reported.}\n"
    @printf moments_fid "\\endfoot\n"

    # Write parameter moments
    sorted_parameters = sort(m.parameters, by = (x -> x.key))
    for param in sorted_parameters
        index = m.keys[param.key]
        @printf moments_fid "\$\%4.99s\$ & " param.tex_label
        @printf moments_fid "%s & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f \\\\\n" outmat[index,:]...
    end

    # Write footer
    write_table_postamble(moments_fid; small=true)

    # Close file
    close(moments_fid)
end

"""
```
prior_posterior_means_table(m, prior_means, post_means; subset_string = "")
```

Produces a table of prior means and posterior means. Saves to:

```
tablespath(m, \"estimate\", \"prior_posterior_means[_sub=\$subset_string].tex\")
```
"""
function prior_posterior_means_table(m::AbstractModel,
                                     prior_means::Vector, post_means::Vector;
                                     subset_string::AbstractString = "")
    # Collect inputs into one matrix
    outmat = [prior_means post_means]

    # Open the TeX file
    basename = "prior_posterior_means"
    if !isempty(subset_string)
        basename *= "_sub=$(subset_string)"
    end
    means_out = tablespath(m, "estimate", "$basename.tex")
    means_fid = open(means_out, "w")

    # Write header
    write_table_preamble(means_fid)
    @printf means_fid "\\vspace*{.5cm}\n"
    @printf means_fid "\\begin{longtable}{ccc}\n"
    @printf means_fid "\\caption{Parameter Estimates: Prior and Posterior Means}\n"
    @printf means_fid "\\\\ \\hline\n"
    @printf means_fid "Parameter & Prior & Posterior\n"
    @printf means_fid "\\\\ \\hline\n"
    @printf means_fid "\\endhead\n"
    @printf means_fid "\\hline\n"
    @printf means_fid "\\endfoot\n"

    # Write results
    sorted_parameters = sort(m.parameters, by = (x -> x.key))
    for param in sorted_parameters
        index = m.keys[param.key]
        @printf means_fid "\$\%4.99s\$ & " param.tex_label
        @printf means_fid "\%8.3f & \%8.3f \\\\\n" outmat[index,:]...
    end

    # Write footer
    write_table_postamble(means_fid)

    # Close file
    close(means_fid)
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
function find_density_bands{T<:AbstractFloat}(draws::Matrix, percent::T; minimize::Bool = true)

    if !(0 < percent <= 1)
        error("percent must be between 0 and 1")
    end

    ndraws, nperiods = size(draws)

    if ndraws == 1
        band = repmat(draws, 2, 1)
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

        band[2,i] = draw_variable_i[low]
        band[1,i] = draw_variable_i[high]
    end

    return band
end

"""
```
find_density_bands{T<:AbstractFloat}(draws::Matrix, percents::Vector{T}; minimize::Bool=true)
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
function find_density_bands{T<:AbstractFloat}(draws::Matrix, percents::Vector{T}; minimize::Bool = true)

    bands = DataFrame()

    for p in percents
        out = find_density_bands(draws, p, minimize = minimize)

        bands[symbol("$(100*p)\% UB")] = vec(out[2,:])
        bands[symbol("$(100*p)\% LB")] = vec(out[1,:])
    end

    bands
end

function write_table_preamble(fid::IOStream)
    @printf fid "\\documentclass[12pt]{article}\n"
    @printf fid "\\usepackage{booktabs}\n"
    @printf fid "\\usepackage[justification=centering]{caption}\n"
    @printf fid "\\usepackage[margin=1in]{geometry}\n"
    @printf fid "\\usepackage{longtable}\n"
    @printf fid "\\begin{document}\n"
    @printf fid "\\pagestyle{empty}\n"
end

# `small`: Whether to print an additional curly bracket after "\end{longtable}" (necessary if
# the table is enclosed by "\small{}")
function write_table_postamble(fid::IOStream; small::Bool=false)
    if small
        @printf fid "\\end{longtable}}\n"
    else
        @printf fid "\\end{longtable}\n"
    end

    @printf fid "\\end{document}"
end