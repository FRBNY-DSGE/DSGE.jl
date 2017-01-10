"""
```
compute_moments(m::AbstractModel, percent::Float64 = 0.90; verbose::Symbol=:none)
```

Computes prior and posterior parameter moments. Tabulates prior mean, posterior mean, and
bands in various LaTeX tables stored `tablespath(m)`.

### Arguments
- `m`: the model object
- `percent`: the percentage of the mass of draws from Metropolis-Hastings included between
  the bands displayed in output tables.
"""
function compute_moments(m::AbstractModel, percent::AbstractFloat = 0.90; verbose::Symbol=:none)

    # Read in the matrix of parameter draws from metropolis-hastings
    filename = rawpath(m, "estimate", "mhsave.h5")

    if VERBOSITY[verbose] >= VERBOSITY[:low]
        @printf "Reading parameter draws from %s\n" filename
    end

    param_draws = []
    try
        h5open(filename, "r") do f
            param_draws = read(f, "mhparams")
        end
    catch
        @printf "Could not open file %s\n" filename
        return
    end

    n_draws = size(param_draws,1)

    # Save mean parameter vector
    save_mean_parameters(m, param_draws)

    # Produce TeX table of moments
    make_moment_tables(m, param_draws, percent, verbose=verbose)
end

"""
```
save_mean_parameters{T<:AbstractFloat}(m::AbstractModel, draws::Matrix{T})
```

Computes and saves the posterior mean of the parameters.

### Arguments
- `m`
- `draws`: n_draws x n_parameters matrix holding the posterior draws from
  Metropolis-Hastings
"""
function save_mean_parameters{T<:AbstractFloat}(m::AbstractModel, draws::Matrix{T})
    post_means = mean(draws,1)'
    filename = workpath(m, "estimate", "paramsmean.h5")
    h5open(filename, "w") do f
        f["params"] = post_means
    end
end

"""
```
make_moment_tables{T<:AbstractFloat}(m::AbstractModel, draws::Matrix{T},
    percent::AbstractFloat; verbose::Symbol = :none)
```

Tabulates parameter moments in 3 LaTeX tables:

1. For MAIN parameters, a list of prior means, prior standard deviations, posterior means,
   and 90% bands for posterior draws
2. For LESS IMPORTANT parameters, a list of the prior means, prior standard deviations,
   posterior means and 90% bands for posterior draws.
3. A list of prior means and posterior means

### Arguments
- `draws`: [n_draws x n_parameters] matrix holding the posterior draws from
  Metropolis-Hastings from save/mhsave.h5
- `percent`: the mass of observations we want; 0 <= percent <= 1
"""
function make_moment_tables{T<:AbstractFloat}(m::AbstractModel, draws::Matrix{T},
                                              percent::AbstractFloat;
                                              verbose::Symbol = :none)

    # STEP 1: Extract moments of prior distribution from m.parameters
    n_params    = length(m.parameters)
    prior_means = zeros(n_params,1)
    prior_std   = zeros(n_params,1)

    # Track prior distributions
    prior_dist  = Vector{AbstractString}(n_params)
    distid(::Distributions.Beta)    = "B"
    distid(::Distributions.Gamma)   = "G"
    distid(::Distributions.Normal)  = "N"
    distid(::DSGE.RootInverseGamma) = "IG"

    for (i,k) in enumerate(m.keys)

        if i > n_params
            continue
        end

        param = getindex(m,i)

        # If param is fixed, we simply set mean to its value and std to 0
        if param.fixed
            prior_means[i]  = param.value
            prior_std[i] = 0
            prior_dist[i] = "-"
        else
            pri = get(param.prior)
            α, β = moments(pri)
            prior_means[i] = α
            prior_std[i] = β
            prior_dist[i] = distid(pri)
        end
    end

    # STEP 2: Compute moments and `percent' bands from parameter draws

    # Posterior mean for each
    post_means = mean(draws,1)'

    # Covariance: Note that this has already been computed and saved
    # in $workpath(m)/parameter_covariance.h5.
    # parameter θ_sig = cov(θ, mean=θ_hat')

    # Bands for each
    post_bands = []
    try
        post_bands = find_density_bands(draws, percent, minimize=true)
    catch
        println("percent must be between 0 and 1")
        return -1
    end

    # We need the transpose
    post_bands = post_bands'

    # STEP 3: Create variables for writing to output table

    # Two-row column names. First row is multicolumn, where entries (i,str) are `i` columns
    # with content `str`.
    colnames0 = [(1,""), (3,"Prior"),(3,"Posterior")]
    colnames = ["Parameter", "Type", "Mean", "SD", "Mean",
                "$(100*percent)\\% {\\tiny Lower Band}",
                "$(100*percent)\\% {\\tiny Upper Band}"]

    # prior type, mean, std dev; posterior mean, 90% lower band, 90% upper bands
    outmat = [prior_dist prior_means prior_std post_means post_bands]

    # prior mean and posterior mean (n_params x 2)
    outmat2 = [prior_means post_means]

    # STEP 4: WRITE TABLES

    # 4a. Write to Table 1: prior mean, std dev, posterior mean and bands

    # Open and start the TeX file
    moments_out = tablespath(m, "estimate", "moments.tex")
    moments_fid = open(moments_out, "w")

    write_table_preamble(moments_fid)

    @printf moments_fid "\\vspace*{.5cm}\n"
    @printf moments_fid "{\\small\n"
    @printf moments_fid "\\begin{longtable}{lcccccc}\n"
    @printf moments_fid "\\caption{Parameter Estimates}\n"
    @printf moments_fid "\\\\ \\hline\n"

    # Column names
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

    # Close the file
    write_table_postamble(moments_fid; small=true)

    # 4b. Write to Table 5: Prior mean and posterior mean

    # Open the TeX file and set up the heading
    means_out = tablespath(m,"estimate", "prior_posterior_means.tex")
    means_fid = open(means_out,"w")

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

    # Write out the results
    for param in sorted_parameters
        index = m.keys[param.key]
        @printf means_fid "\$\%4.99s\$ & " param.tex_label
        @printf means_fid "\%8.3f & \%8.3f \\\\\n" outmat2[index,:]...
    end

    write_table_postamble(means_fid)

    if VERBOSITY[verbose] >= VERBOSITY[:low]
        @printf "Tables are saved as %s.\n" tablespath(m, "estimate", "*.tex")
    end
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
