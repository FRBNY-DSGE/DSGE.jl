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

    # 4a. Write to Table 1: prior mean, std dev, posterior mean and bands for IMPORTANT
    #     parameters

    # Open and start the TeX file
    moments0_out = tablespath(m, "estimate", "moments0.tex")
    moments0_fid = open(moments0_out, "w")

    write_table_preamble(moments0_fid)

    @printf moments0_fid "\\caption{Parameter Estimates}\n"
    @printf moments0_fid "\\vspace*{.5cm}\n"
    @printf moments0_fid "{\\small \n"
    @printf moments0_fid "\\begin{tabular}{lcccccc}\n"
    @printf moments0_fid "\\hline\n"

    # Column names
    @printf moments0_fid "\\multicolumn{%d}{c}{%s}" colnames0[1][1] colnames0[1][2]
    for col in colnames0[2:end]
        @printf moments0_fid " & \\multicolumn{%d}{c}{%s}" col[1] col[2]
    end
    @printf moments0_fid "\\\\ \n"
    @printf moments0_fid "%s" colnames[1]
    for col in colnames[2:end]
        @printf moments0_fid " & %s" col
    end
    @printf moments0_fid "\\\\ \n"
    @printf moments0_fid "\\cmidrule(lr){1-1}\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}\n"

    # Keep track of indices for important parameters
    important_para = []

    for (index, param) in enumerate(m.parameters)

        if (!ismatch(r"rho_", param.tex_label) &&
            !ismatch(r"zeta_", param.tex_label) &&
            !ismatch(r"psi_", param.tex_label) &&
            !ismatch(r"nu_l", param.tex_label) &&
            !ismatch(r"pi\^\*", param.tex_label) &&
            !ismatch(r"sigma_{pi}\^\*",param.tex_label) &&
            (!ismatch(r"pistar", param.tex_label)))
            continue
        end

        if ismatch(r"rho_chi",param.tex_label)
            continue
        end

        #Print the parameter name and values in outmat
        @printf moments0_fid "\$\%4.99s\$ & " param.tex_label
        @printf moments0_fid "%s & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f \\\\\n" outmat[index,:]...

        important_para = [important_para; index]
    end

    # Close the file
    write_table_postamble(moments0_fid;
        note="{\\footnotesize Note: For Inverse Gamma (IG) prior mean and SD, \$\\tau^2\$ and \$\\nu\$
        reported.}",
        small=true)

    # 4b. Write to Table 2: Prior mean, std dev and posterior mean, bands for other params
    moments_table_out = tablespath(m, "estimate", "moments1.tex")
    moments_table_fid = open(moments_table_out, "w")

    write_table_preamble(moments_table_fid)

    @printf moments_table_fid "\\caption{Parameter Estimates}\n"
    @printf moments_table_fid "\\vspace*{.2cm}\n"
    @printf moments_table_fid "{\\small \n"
    @printf moments_table_fid "\\begin{tabular}{lcccccc}\n"
    @printf moments_table_fid "\\hline\n"

    # Column names
    for col in colnames0
        @printf moments_table_fid "\\multicolumn{%d}{c}{%s} &" col[1] col[2]
    end
    @printf moments_table_fid "\\\n"
    for col in colnames
        @printf moments_table_fid "%s & " col
    end
    @printf moments_table_fid "\\\n"
    @printf moments_table_fid "\\cmidrule(lr){1-1}\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}\n"

    # Counter for parameters to track length of table and number of tables in excess of
    # default (1)
    other_para  = 1
    table_count = 0

    for (index, param) in enumerate(m.parameters)

        if in(index, important_para)
            continue
        end

        # Make a new table if the current one is too large
        if other_para % 25 == 0 && index ≠ length(m.parameters)

            # Close the file
            write_table_postamble(moments_table_fid;
                note="{\\footnotesize Note: For Inverse Gamma (IG) prior mean and SD, \$\\tau^2\$ and \$\\nu\$
                reported.}",
                small=true)

            # Update table counter
            table_count += 1

            # Start the new file
            filename = @sprintf "moments%d.tex" table_count
            moments_table_out = tablespath(m,"estimate",filename)
            moments_table_fid = open(moments_table_out,"w")

            write_table_preamble(moments_table_fid)
            @printf moments_table_fid "\\caption{Parameter Estimates}\n"
            @printf moments_table_fid "\\vspace*{.2cm}\n"
            @printf moments_table_fid "{\\small \n"
            @printf moments_table_fid "\\begin{tabular}{lcccccc}\n"
            @printf moments_table_fid "\\hline\n"

            # Column names
            @printf moments_table_fid "\\multicolumn{%d}{c}{%s}" colnames0[1][1] colnames0[1][2]
            for col in colnames0[2:end]
                @printf moments_table_fid " & \\multicolumn{%d}{c}{%s}" col[1] col[2]
            end
            @printf moments_table_fid "\\\\ \n"
            @printf moments_table_fid "%s" colnames[1]
            for col in colnames[2:end]
                @printf moments_table_fid " & %s" col
            end
            @printf moments_table_fid "\\\\ \n"
            @printf moments_table_fid "\\cmidrule(lr){1-1}\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}\n"
        end

        #Print the parameter name and values in outmat
        @printf moments_table_fid "\$\%4.99s\$ & " param.tex_label
        @printf moments_table_fid "%s & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f \\\\\n" outmat[index,:]...

        other_para += 1
    end

    write_table_postamble(moments_table_fid; small=true)

    # 4c. Write to Table 5: Prior mean and posterior mean for all parameters

    # Keep track of how many tables we've made
    table_count = 0

    # Open the TeX file and set up the heading
    prioPostMean_out = tablespath(m,"estimate", "moments_prioPostMean.tex")
    prioPostMean_fid = open(prioPostMean_out,"w")

    write_table_preamble(prioPostMean_fid)
    @printf prioPostMean_fid "\\caption{Parameter Estimates: Prior and Posterior Mean}\n"
    @printf prioPostMean_fid "\\vspace*{.5cm}\n"
    @printf prioPostMean_fid "\\begin{tabular}{ccc}\\hline \n"
    @printf prioPostMean_fid " Parameter & Prior & Posterior  \\tabularnewline \\hline\n"

    # Write out the results
    for (index, param) in enumerate(m.parameters)

        if index % 40 == 0 && index ≠ length(m.parameters)

            # Close the old file
            write_table_postamble(prioPostMean_fid)

            # Generate the new filename
            table_count += 1
            filename = @sprintf "moments_prioPostMean_%d.tex" table_count

            # Open a new file and start the next table
            prioPostMean_out = tablespath(m,"estimate",filename)
            prioPostMean_fid = open(prioPostMean_out,"w")
            write_table_preamble(prioPostMean_fid)
            @printf prioPostMean_fid "\\caption{Parameter Estimates: Prior and Posterior Mean}\n"
            @printf prioPostMean_fid "\\vspace*{.5cm}\n"
            @printf prioPostMean_fid "\\begin{tabular}{ccc}\\hline \n"
            @printf prioPostMean_fid " Parameter & Prior & Posterior  \\tabularnewline \\hline\n"
        end

        @printf prioPostMean_fid "\$\%4.99s\$ & " param.tex_label

        val = outmat2[index,:]
        @printf prioPostMean_fid "\%8.3f & \%8.3f \\\\\n" val[1] val[2]
    end

    write_table_postamble(prioPostMean_fid)

    if VERBOSITY[verbose] >= VERBOSITY[:low]
        @printf "Tables are saved as %s.\n" tablespath(m, "estimate", "*.tex")
    end
end

function write_table_preamble(fid::IOStream)
    @printf fid "\\documentclass[12pt]{article}\n"
    @printf fid "\\usepackage{booktabs}\n"
    @printf fid "\\begin{document}\n"
    @printf fid "\\pagestyle{empty}\n"
    @printf fid "\\begin{table}[h]\n"
    @printf fid "\\centering\n"
end


# `small`: Whether to print an additional curly bracket after "\end{tabular}" (necessary if
# the table is enclosed by "\small{}")

function write_table_postamble(fid::IOStream; note::AbstractString="", small::Bool=false)
    @printf fid "\\\\ \\\hline\n"

    if small
        @printf fid "\\end{tabular}}\n"
    else
        @printf fid "\\end{tabular}\n"
    end

    if !isempty(note)
        @printf fid "%s\n" note
    end

    @printf fid "\\end{table}\n"
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
- `draws`: Matrix of parameter draws (from Metropolis-Hastings, for example)
- `percent`: percent of data within bands (e.g. .9 to get 90% of mass within bands)

### Optional Arguments
- `minimize`: if `true`, choose shortest interval, otherwise just chop off lowest and
  highest (percent/2)
"""
function find_density_bands(draws::Matrix, percent::AbstractFloat; minimize::Bool=true)

    if !(0 < percent < 1)
        error("percent must be between 0 and 1")
    end

    n_draws, n_draw_dimensions = size(draws)
    band = zeros(2, n_draw_dimensions)
    n_in_band  = round(Int, percent * n_draws)

    for i in 1:n_draw_dimensions

        # Sort response for parameter i such that 1st element is largest
        draw_variable_i = draws[:,i]
        sort!(draw_variable_i, rev=true)

        # Search to find the interval containing the minimum # of observations
        # comprising `percent` of the mass
        if minimize

            upper_index = 1
            done        = 0
            j           = 2
            minwidth = draw_variable_i[1] - draw_variable_i[n_in_band]

            while j <= (n_draws - n_in_band + 1)

                newwidth = draw_variable_i[j] - draw_variable_i[j + n_in_band - 1]

                if newwidth < minwidth
                    upper_index = j
                    minwidth = newwidth
                end

                j += 1
            end

        else
            upper_index = n_draws - nwidth - floor(.5*n_draws-n_in_band)
        end

        band[2,i] = draw_variable_i[upper_index]
        band[1,i] = draw_variable_i[upper_index + n_in_band - 1]
    end

    return band
end
