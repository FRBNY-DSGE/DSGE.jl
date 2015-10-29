## moments.jl: Computes and tabulates moments of parameter draws from Metropolis-Hastings

using HDF5, Compat

#=
doc"""
compute_moments{T<:AbstractDSGEModel}(m::T, percent::Float64 = 0.90)

### Parameters
  - `m`: the model object
  - `percent`: the percentage of the mass of draws from Metropolis-Hastings included between the bands displayed in output tables. 

### Description
Computes prior and posterior parameter moments. Tabulates prior mean, posterior mean, and
bands in various LaTeX tables stored `tablespath(m)`.
"""
=#
function compute_moments{T<:AbstractDSGEModel}(m::T, percent::Float64 = 0.90; verbose=true)
    
    ### Step 1: Read in the matrix of parameter draws from metropolis-hastings

    filename = joinpath(rawpath(m,"estimate"),"sim_save.h5")

    if verbose
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

    num_draws = size(param_draws,1)

    
    ### Step 2: Produce TeX table of moments

    make_moment_tables(m,param_draws,percent, verbose=verbose)

end

#=
doc"""
make_moment_tables{T<:AbstractFloat}(m::AbstractDSGEModel, Θ::Array{T,2}, percent::Float64)

### Parameters
    - `Θ`: [num_draws x num_parameters] matrix holding the posterior draws from metropolis-hastings
           from save/sim_save.h5 
    - `percent`: the mass of observations we want; 0 <= percent <= 1

### Description
Tabulates parameter moments in 3 LaTeX tables:
    
1. For MAIN parameters, a list of prior means, prior standard deviations, posterior means, and 90% bands for posterior draws

2. For LESS IMPORTANT parameters, a list of the prior means, prior standard deviations, posterior means and 90% bands for posterior draws.

3. A list of prior means and posterior means
"""
=#
function make_moment_tables{T<:AbstractFloat}(m::AbstractDSGEModel, Θ::Array{T,2}, percent::Float64; verbose=true)
    

    ###########################################################################################
    ## STEP 1: Extract moments of prior distribution from m.parameters
    ###########################################################################################
    
    num_params = length(m.parameters) 
    prior_means = zeros(num_params,1)
    prior_stddev = zeros(num_params,1)

    
    for (i,k) in enumerate(m.keys)

        if(i > num_params)
            continue
        end
        
        param  = getindex(m,i)
        
        if isa(param.prior.value, DSGE.Normal)

            prior_means[i] = param.prior.value.μ
            prior_stddev[i] = param.prior.value.σ
            
        elseif isa(param.prior.value, Distributions.Beta)
            μ,σ = betaMoments(param.prior.value)
            
            prior_means[i] = μ
            prior_stddev[i] = σ

        elseif isa(param.prior.value, Distributions.Gamma)
            μ,σ = gammaMoments(param.prior.value)
            
            prior_means[i] = μ
            prior_stddev[i] = σ  # small \theta
            
        end
    end
    
    ###########################################################################################
    ## STEP 2: Compute moments and `percent' bands from parameter draws
    ###########################################################################################

    # Posterior mean for each
    Θ_hat = mean(Θ,1)'       

    # Covariance: Note that this has already been computed and saved
    # in $outpath(m)/parameter_covariance.h5.   
    # parameter Θ_sig = cov(Θ, mean=Θ_hat') 

    # Bands for each
    Θ_bands = []

    try
        Θ_bands = find_density_bands(Θ,percent,minimize=true)' # We need the transpose
    catch
        println("percent must be between 0 and 1")
        return -1
    end
    

    ###########################################################################################
    ## STEP 3: Create variables for writing to output table
    ###########################################################################################
    
    colnames = ["Parameter ", "Prior Mean ", "Prior Stdd ", "Post Mean ",
                "$(100*percent)\\% {\\tiny Lower Band} ", "$(100*percent)\\% {\\tiny Upper Band} " ]

    outmat = [prior_means prior_stddev Θ_hat Θ_bands]      # prior mean and std dev, posterior mean,
                                                           # and bands (n_params x 5)

    outmat2 = [prior_means Θ_hat]                          # prior mean and posterior mean
                                                           # (n_params x 2)


    ###########################################################################################
    ## STEP 4: WRITE TABLES
    ###########################################################################################
    
    ###########################################################################################
    ## 4a. Write to Table 1: prior mean, std dev, posterior mean and bands for IMPORTANT parameters
    ###########################################################################################
    
    # Open and start the TeX file
    mainParams_out = joinpath(tablespath(m,"estimate"), "moments_mainParams.tex")
    mainParams_fid = open(mainParams_out,"w")

    beginTexTableDoc(mainParams_fid)

    @printf(mainParams_fid,"\\caption{Parameter Estimates}\n")
    @printf(mainParams_fid,"\\vspace*{.5cm}\n")
    @printf(mainParams_fid,"{\\small \n")
    @printf(mainParams_fid,"\\begin{tabular}{lllllll}\\hline \n")

    # Column names
    for col in colnames
        @printf(mainParams_fid, "%4.99s & ", col)
    end
        
    # Keep track of indices for important parameters 
    important_para = []

    for (index, param) in enumerate(m.parameters)   
        
        if (!ismatch(r"rho_", param.texLabel) &&
            !ismatch(r"zeta_", param.texLabel) &&
            !ismatch(r"psi_", param.texLabel) &&
            !ismatch(r"nu_l", param.texLabel) &&
            !ismatch(r"pi\^\*", param.texLabel) &&
            !ismatch(r"sigma_{pi}\^\*",param.texLabel) &&
            (!ismatch(r"pistar", param.texLabel)))
            # (!ismatch(r"ups", param.texLabel)))  ##Is this correct? mspec == 16 or u_^*
            continue
        end
            
        # TODO: Decide whether subspec should be a field in the model
        if(ismatch(r"rho_chi",param.texLabel)) # ??? || (isequal(subspec,7) && texLabel == ":rho_b"))
            continue
        end

        @printf(mainParams_fid, "\\\\ \n \$\%4.99s\$ & ", param.texLabel)
        
        #Print the values in outmat
        for val in outmat[index,:]
            @printf(mainParams_fid, "\%8.3f & ",val)
        end

        important_para = [important_para; index]
        
    end

    # Close the file
    endTexTableDoc(mainParams_fid;small=true)

    ###########################################################################################
    ## 4b. Write to Table 2: Prior mean, std dev and posterior mean, bands for other params
    ###########################################################################################
    
    periphParams_out = joinpath(tablespath(m,"estimate"), "moments_periphParams_0.tex")
    periphParams_fid = open(periphParams_out,"w")

    beginTexTableDoc(periphParams_fid)

    @printf(periphParams_fid,"\\caption{Parameter Estimates}\n")
    @printf(periphParams_fid,"\\vspace*{.2cm}\n")
    @printf(periphParams_fid,"{\\small \n")
    @printf(periphParams_fid,"\\begin{tabular}{lllllll}\\hline \n")
    
    # Column names
    for col in colnames
        @printf(periphParams_fid, "%4.99s & ", col)
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
            endTexTableDoc(periphParams_fid;small=true)

            # Update table counter
            table_count = table_count + 1

            # Start the new file
            filename = @sprintf("moments_periphParams_%d.tex",table_count)
            periphParams_out = joinpath(tablespath(m,"estimate"),filename)
            periphParams_fid = open(periphParams_out,"w")
            
            beginTexTableDoc(periphParams_fid)
            @printf(periphParams_fid,"\\caption{Parameter Estimates}\n")
            @printf(periphParams_fid,"\\vspace*{.2cm}\n")
            @printf(periphParams_fid,"{\\small \n")
            @printf(periphParams_fid,"\\begin{tabular}{lllllll}\\hline \n")
            
            for col in colnames
                @printf(periphParams_fid, "%4.99s & ", col)
            end
        end
        
        # Print the parameter name
        @printf(periphParams_fid, "\\\\ \n \$\%4.99s\$ & ", param.texLabel)
        
        # Print the values in outmat
        for val in outmat[index,:]
            @printf(periphParams_fid, "\%8.3f & ",val)
        end
        
        other_para = other_para+1
        
    end
    endTexTableDoc(periphParams_fid;small=true)

    ###########################################################################################
    ## 4c. Write to Table 5: Prior mean and posterior mean for all parameters
    ###########################################################################################

    table_count = 0  # Keep track of how many tables we've made
    
    # Open the TeX file and set up the heading
    prioPostMean_out = joinpath(tablespath(m,"estimate"), "moments_prioPostMean.tex")
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
            prioPostMean_out = joinpath(tablespath(m,"estimate"),filename)
            prioPostMean_fid = open(prioPostMean_out,"w")            
            beginTexTableDoc(prioPostMean_fid)
            @printf(prioPostMean_fid,"\\caption{Parameter Estimates: Prior and Posterior Mean}\n")
            @printf(prioPostMean_fid,"\\vspace*{.5cm}\n")
            @printf(prioPostMean_fid,"\\begin{tabular}{ccc}\\hline \n")
            @printf(prioPostMean_fid," Parameter & Prior & Posterior  \\tabularnewline \\hline\n")

        end

        
        @printf(prioPostMean_fid, "\\\\ \n \$\%4.99s\$ & ", param.texLabel)

        val = outmat2[index,:]
        @printf(prioPostMean_fid, "\%8.3f &   \%8.3f  ", val[1], val[2])
        
    end
    
    endTexTableDoc(prioPostMean_fid)

    if verbose
        println("Tables are in ",tablespath(m,"estimate"))
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

#=
doc"""
### Parameters
- `fid`: File descriptor

### Optional Arguments
- `small`: Whether to print an additional curly bracket after "\end{tabular}" (necessary if the table is enclosed by "\small{}")

### Description
Prints the necessarily lines to end a table and close a LaTeX document to file descriptor `fid`, then closes the file.
"""
=#
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

#=
doc"""
find_density_bands(draws::Matrix, percent::Real; minimize::Bool=true)

### Parameters
- draws: Matrix of parameter draws (from Metropolis-Hastings, for example)
- percent: percent of data within bands (e.g. .9 to get 90% of mass within bands)

### Optional Arguments
- `minimize`: if `true`, choose shortest interval, otherwise just chop off lowest and highest (percent/2)

### Description
Returns a [2 x cols(draws)] matrix `bands` such that `percent` of the mass of `draws[:,i]` is above `bands[1,i]` and below `bands[2,i]`.
"""
=#
function find_density_bands(draws::Matrix, percent::Real; minimize::Bool=true)

    if(percent < 0 || percent > 1)
        error("percent must be between 0 and 1")
    end
    
    num_draws, num_draw_dimensions = size(draws)
    band  = zeros(2, num_draw_dimensions)
    num_in_band  = round(Int, percent * num_draws)
    
    for i in 1:num_draw_dimensions

        # Sort response for parameter i such that 1st element is largest
        draw_variable_i = draws[:,i]
        sort!(draw_variable_i, rev=true)

        # Search to find the interval containing the minimum # of observations
        # comprising `percent` of the mass
        if minimize

            upper_index=1
            minwidth = draw_variable_i[1] - draw_variable_i[num_in_band]
            done = 0
            j = 2
            
            while j <= (num_draws - num_in_band + 1)

                newwidth = draw_variable_i[j] - draw_variable_i[j + num_in_band - 1]

                if newwidth < minwidth
                    upper_index = j
                    minwidth = newwidth        
                end

                j = j+1
            end
            
        else
            upper_index = num_draws - nwidth - floor(.5*num_draws-num_in_band)
        end

        band[2,i] = draw_variable_i[upper_index]
        band[1,i] = draw_variable_i[upper_index + num_in_band - 1]
    end

    return band
end
