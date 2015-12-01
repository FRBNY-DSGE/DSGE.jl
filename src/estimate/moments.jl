"""
`compute_moments{T<:AbstractModel}(m::T, percent::Float64 = 0.90; verbose::Symbol=:none)`

Computes prior and posterior parameter moments. Tabulates prior mean, posterior mean, and
bands in various LaTeX tables stored `tablespath(m)`.

### Arguments
  - `m`: the model object
  - `percent`: the percentage of the mass of draws from Metropolis-Hastings included between
    the bands displayed in output tables. 
"""
function compute_moments{T<:AbstractModel}(m::T, percent::Float64 = 0.90; 
                                               verbose::Symbol=:none)
    
    ### Step 1: Read in the matrix of parameter draws from metropolis-hastings

    filename = rawpath(m,"estimate","mh_save.h5")

    if VERBOSITY[verbose] >= VERBOSITY[:low]
        println("Reading draws from Metropolis-Hastings from $filename...")
    end
    
    param_draws = []
    post = []

    try
        fid = h5open(filename,"r+") 
        param_draws = read(fid,"parasim")
        #post = read(fid,"postsim")
        close(fid)
    catch
        @printf(1,"Could not open file %s", filename)
    end

    n_draws = size(param_draws,1)

    
    ### Step 2: Produce TeX table of moments

    make_moment_tables(m,param_draws,percent, verbose=verbose)

end

"""
```
make_moment_tables{T<:AbstractFloat}(m::AbstractModel, θ::Array{T,2}, percent::Float64;
                                     verbose::Symbol=:none)
```

Tabulates parameter moments in 3 LaTeX tables:
    
1. For MAIN parameters, a list of prior means, prior standard deviations, posterior means,
   and 90% bands for posterior draws

2. For LESS IMPORTANT parameters, a list of the prior means, prior standard deviations,
   posterior means and 90% bands for posterior draws.

3. A list of prior means and posterior means

### Arguments
    - `θ`: [n_draws x n_parameters] matrix holding the posterior draws from metropolis-hastings
           from save/mh_save.h5 
    - `percent`: the mass of observations we want; 0 <= percent <= 1
"""
function make_moment_tables{T<:AbstractFloat}(m::AbstractModel,
                                              θ::Array{T,2},
                                              percent::Float64;
                                              verbose::Symbol=:none)

    ########################################################################################
    ## STEP 1: Extract moments of prior distribution from m.parameters
    ########################################################################################
    
    n_params = length(m.parameters) 
    prior_means = zeros(n_params,1)
    prior_std = zeros(n_params,1)

    
    for (i,k) in enumerate(m.keys)

        if(i > n_params)
            continue
        end
        
        param  = getindex(m,i)
        
        # Null prior means parameter is fixed
        if isnull(param.prior)
            prior_means[i]  = param.value
            prior_std[i] = 0

        # TODO xxxMoments could be declared something like
        # function distMoments{T<:Distribution}(prior_value::T)
        elseif isa(param.prior.value, DSGE.Normal)

            prior_means[i] = param.prior.value.μ
            prior_std[i] = param.prior.value.σ
            
        elseif isa(param.prior.value, Distributions.Beta)
            μ,σ = moments(param.prior.value)
            
            prior_means[i] = μ
            prior_std[i] = σ

        elseif isa(param.prior.value, Distributions.Gamma)
            μ,σ = moments(param.prior.value)
            
            prior_means[i] = μ
            prior_std[i] = σ  # small \theta
            
        end
    end
    
    ########################################################################################
    ## STEP 2: Compute moments and `percent' bands from parameter draws
    ########################################################################################

    # Posterior mean for each
    θ_hat = mean(θ,1)'       

    # Covariance: Note that this has already been computed and saved
    # in $workpath(m)/parameter_covariance.h5.   
    # parameter θ_sig = cov(θ, mean=θ_hat') 

    # Bands for each
    θ_bands = []

    try
        θ_bands = find_density_bands(θ,percent,minimize=true)' # We need the transpose
    catch
        println("percent must be between 0 and 1")
        return -1
    end
    

    ########################################################################################
    ## STEP 3: Create variables for writing to output table
    ########################################################################################
    
    colnames = ["Parameter ", "Prior Mean ", "Prior Stdd ", "Post Mean ",
                "$(100*percent)\\% {\\tiny Lower Band} ", "$(100*percent)\\% {\\tiny Upper Band} " ]

    outmat = [prior_means prior_std θ_hat θ_bands]      # prior mean and std dev, posterior mean,
                                                           # and bands (n_params x 5)

    outmat2 = [prior_means θ_hat]                          # prior mean and posterior mean
                                                           # (n_params x 2)


    ########################################################################################
    ## STEP 4: WRITE TABLES
    ########################################################################################
    
    ########################################################################################
    ## 4a. Write to Table 1: prior mean, std dev, posterior mean and bands for IMPORTANT
    ##     parameters
    ########################################################################################
    
    # Open and start the TeX file
    moments0_out = tablespath(m,"estimate", "moments0.tex")
    moments0_fid = open(moments0_out,"w")

    beginTexTableDoc(moments0_fid)

    @printf(moments0_fid,"\\caption{Parameter Estimates}\n")
    @printf(moments0_fid,"\\vspace*{.5cm}\n")
    @printf(moments0_fid,"{\\small \n")
    @printf(moments0_fid,"\\begin{tabular}{lllllll}\\hline \n")

    # Column names
    for col in colnames
        @printf(moments0_fid, "%4.99s & ", col)
    end
        
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
            # (!ismatch(r"ups", param.tex_label)))  ##Is this correct? mspec == 16 or u_^*
            continue
        end
            
        # TODO: Decide whether subspec should be a field in the model
        if(ismatch(r"rho_chi",param.tex_label)) # ??? || (isequal(subspec,7) && tex_label == ":rho_b"))
            continue
        end

        @printf(moments0_fid, "\\\\ \n \$\%4.99s\$ & ", param.tex_label)
        
        #Print the values in outmat
        for val in outmat[index,:]
            @printf(moments0_fid, "\%8.3f & ",val)
        end

        important_para = [important_para; index]
        
    end

    # Close the file
    endTexTableDoc(moments0_fid;small=true)

    ########################################################################################
    ## 4b. Write to Table 2: Prior mean, std dev and posterior mean, bands for other params
    ########################################################################################
    
    moments_table_out = tablespath(m,"estimate", "moments1.tex")
    moments_table_fid = open(moments_table_out,"w")

    beginTexTableDoc(moments_table_fid)

    @printf(moments_table_fid,"\\caption{Parameter Estimates}\n")
    @printf(moments_table_fid,"\\vspace*{.2cm}\n")
    @printf(moments_table_fid,"{\\small \n")
    @printf(moments_table_fid,"\\begin{tabular}{lllllll}\\hline \n")
    
    # Column names
    for col in colnames
        @printf(moments_table_fid, "%4.99s & ", col)
    end
    
    # Counter for parameters to track length of table and number of tables in excess of
    # default (1)
    other_para = 1
    table_count = 0

    for (index, param) in enumerate(m.parameters)
    
        if in(index, important_para)
            continue
        end
        
        # Make a new table if the current one is too large
        
        if ((other_para%25 == 0) && (index != length(m.parameters)) )

            # Finish and close the old file
            endTexTableDoc(moments_table_fid;small=true)

            # Update table counter
            table_count = table_count + 1

            # Start the new file
            filename = @sprintf("moments%d.tex",table_count)
            moments_table_out = tablespath(m,"estimate",filename)
            moments_table_fid = open(moments_table_out,"w")
            
            beginTexTableDoc(moments_table_fid)
            @printf(moments_table_fid,"\\caption{Parameter Estimates}\n")
            @printf(moments_table_fid,"\\vspace*{.2cm}\n")
            @printf(moments_table_fid,"{\\small \n")
            @printf(moments_table_fid,"\\begin{tabular}{lllllll}\\hline \n")
            
            for col in colnames
                @printf(moments_table_fid, "%4.99s & ", col)
            end
        end
        
        # Print the parameter name
        @printf(moments_table_fid, "\\\\ \n \$\%4.99s\$ & ", param.tex_label)
        
        # Print the values in outmat
        for val in outmat[index,:]
            @printf(moments_table_fid, "\%8.3f & ",val)
        end
        
        other_para = other_para+1
        
    end
    endTexTableDoc(moments_table_fid;small=true)

    ########################################################################################
    ## 4c. Write to Table 5: Prior mean and posterior mean for all parameters
    ########################################################################################

    table_count = 0  # Keep track of how many tables we've made
    
    # Open the TeX file and set up the heading
    prioPostMean_out = tablespath(m,"estimate", "moments_prioPostMean.tex")
    prioPostMean_fid = open(prioPostMean_out,"w")

    beginTexTableDoc(prioPostMean_fid)
    @printf(prioPostMean_fid,"\\caption{Parameter Estimates: Prior and Posterior Mean}\n")
    @printf(prioPostMean_fid,"\\vspace*{.5cm}\n")
    @printf(prioPostMean_fid,"\\begin{tabular}{ccc}\\hline \n")
    @printf(prioPostMean_fid," Parameter & Prior & Posterior  \\tabularnewline \\hline\n")
    
    # Write out the results
    for (index, param) in enumerate(m.parameters)

        if ((index % 40 == 0) && (index != length(m.parameters)) )

            # Close the old file
            endTexTableDoc(prioPostMean_fid)

            # Generate the new filename
            table_count = table_count + 1
            filename = @sprintf("moments_prioPostMean_%d.tex",table_count)

            # Open a new file and start the next table
            prioPostMean_out = tablespath(m,"estimate",filename)
            prioPostMean_fid = open(prioPostMean_out,"w")            
            beginTexTableDoc(prioPostMean_fid)
            @printf(prioPostMean_fid,"\\caption{Parameter Estimates: Prior and Posterior Mean}\n")
            @printf(prioPostMean_fid,"\\vspace*{.5cm}\n")
            @printf(prioPostMean_fid,"\\begin{tabular}{ccc}\\hline \n")
            @printf(prioPostMean_fid," Parameter & Prior & Posterior  \\tabularnewline \\hline\n")

        end

        
        @printf(prioPostMean_fid, "\\\\ \n \$\%4.99s\$ & ", param.tex_label)

        val = outmat2[index,:]
        @printf(prioPostMean_fid, "\%8.3f &   \%8.3f  ", val[1], val[2])
        
    end
    
    endTexTableDoc(prioPostMean_fid)

    if VERBOSITY[verbose] >= VERBOSITY[:low]
        @printf "Tables are saved as %s.\n" tablespath(m, "estimate", "*.tex")
    end
end

#=
doc"""
beginTexTableDoc(fid::IOStream)

### Parameters
- `fid`: File descriptor 

### Description
Prints the preamble for a LaTeX table to the file indicated by `fid`.
"""
=#
function beginTexTableDoc(fid::IOStream)

    @printf(fid,"\\documentclass[12pt]{article}\n")
    @printf(fid,"\\usepackage[dvips]{color}\n")
    @printf(fid,"\\begin{document}\n")
    @printf(fid,"\\pagestyle{empty}\n")
    @printf(fid,"\\begin{table}[h] \\centering\n")
    
end

"""
`endTexTableDoc(fid::IOStream;small::Bool=false)`

Prints the necessarily lines to end a table and close a LaTeX document to file descriptor `fid`, then closes the file.

### Arguments
- `fid`: File descriptor

### Optional Arguments
- `small`: Whether to print an additional curly bracket after "\end{tabular}" (necessary if the table is enclosed by "\small{}")
"""
function endTexTableDoc(fid::IOStream;small::Bool=false)

    @printf(fid, "\\\\ \\\hline\n")
    
    if small
        @printf(fid,"\\end{tabular}}\n")
    else
        @printf(fid,"\\end{tabular}\n")
    end
    
    @printf(fid,"\\end{table}\n")
    @printf(fid,"\\end{document}")
    close(fid)

end

"""
`find_density_bands(draws::Matrix, percent::Real; minimize::Bool=true)`

### Parameters
- draws: Matrix of parameter draws (from Metropolis-Hastings, for example)
- percent: percent of data within bands (e.g. .9 to get 90% of mass within bands)

### Optional Arguments
- `minimize`: if `true`, choose shortest interval, otherwise just chop off lowest and
  highest (percent/2)

### Description
Returns a [2 x cols(draws)] matrix `bands` such that `percent` of the mass of `draws[:,i]`
is above `bands[1,i]` and below `bands[2,i]`.
"""
function find_density_bands(draws::Matrix, percent::Real; minimize::Bool=true)

    if(percent < 0 || percent > 1)
        error("percent must be between 0 and 1")
    end
    
    n_draws, n_draw_dimensions = size(draws)
    band  = zeros(2, n_draw_dimensions)
    n_in_band  = round(Int, percent * n_draws)
    
    for i in 1:n_draw_dimensions

        # Sort response for parameter i such that 1st element is largest
        draw_variable_i = draws[:,i]
        sort!(draw_variable_i, rev=true)

        # Search to find the interval containing the minimum # of observations
        # comprising `percent` of the mass
        if minimize

            upper_index=1
            minwidth = draw_variable_i[1] - draw_variable_i[n_in_band]
            done = 0
            j = 2
            
            while j <= (n_draws - n_in_band + 1)

                newwidth = draw_variable_i[j] - draw_variable_i[j + n_in_band - 1]

                if newwidth < minwidth
                    upper_index = j
                    minwidth = newwidth        
                end

                j = j+1
            end
            
        else
            upper_index = n_draws - nwidth - floor(.5*n_draws-n_in_band)
        end

        band[2,i] = draw_variable_i[upper_index]
        band[1,i] = draw_variable_i[upper_index + n_in_band - 1]
    end

    return band
end
