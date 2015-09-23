## moments.jl: Computes and tabulates moments of parameter draws from Metropolis-Hastings

using HDF5, Compat

function computeMoments{T<:AbstractDSGEModel}(m::T, percent::Float64 = 0.90)
    
    ## Computes prior and posterior parameter moments, tabulates them in various TeX tables,
    ## and plots parameter draws from the prior and posterior distribution
    ##
    ## Inputs:
    ## - m: an instance of an AbstractDSGEModel subtype
    ## - percent: the percentage of the mass of draws from Metropolis-Hastings
    ##            we want to include between bands shown in the table

    
    # Read in the matrix of parameter draws from metropolis-hastings

    infile = joinpath(outpath(m),"sim_save.h5")
    println("Reading draws from Metropolis-Hastings from $infile...")

    Θ = []
    post = []

    try
        fid = h5open(infile,"r") do fid
            Θ = read(fid,"parasim")
            #post = read(fid,"postsim")
        end
    catch
        @printf(1,"Could not open file %s", infile)
    end

    # Convert back to Float64 for compatability with other variables
    #Θ = convert(Matrix{Float64},Θ)

    num_draws = size(Θ,1)

    # Produce TeX table of moments
    makeMomentTables(m,Θ,percent)

end


function makeMomentTables{T<:AbstractFloat}(m::AbstractDSGEModel, Θ::Array{T,2}, percent::Float64)
    
    ## Tabulates parameter moments in 3 LaTeX tables:
    ##
    ## -For MAIN parameters, a list of prior means, prior standard deviations, posterior means,
    ## and 90% bands for posterior draws
    ## -For LESS IMPORTANT parameters, a list of the prior means, prior standard deviations,
    ## posterior means and 90% bands for posterior draws.
    ## -A list of prior means and posterior means
    ##
    ## Input:
    ## - `Θ`: [num_draws x num_parameters] matrix holding the posterior draws from metropolis-hastings
    ##        from save/sim_save.h5 
    ## - `percent`: the mass of observations we want; 0 <= percent <= 1


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
        
        if isa(param.priordist, DSGE.Normal)

            prior_means[i] = param.priordist.μ
            prior_stddev[i] = param.priordist.σ
            
        elseif isa(param.priordist, Distributions.Beta)
            μ,σ = betaMoments(param.priordist)
            
            prior_means[i] = μ
            prior_stddev[i] = σ

        elseif isa(param.priordist, Distributions.Gamma)
            μ,σ = gammaMoments(param.priordist)
            
            prior_means[i] = μ
            prior_stddev[i] = σ  # small \theta
            
        end
    end
    
    ###########################################################################################
    ## STEP 2: Compute moments and `percent' bands from parameter draws
    ###########################################################################################
    
    Θ_hat = mean(Θ,1)'                    # posterior mean for each parameter
    Θ_sig = cov(Θ, mean=Θ_hat')           # posterior std dev for each parameter
    Θ_bands = []
    
    # TODO: Do we want to be doing error handling like this?
    try
        Θ_bands = find_density_bands(Θ,percent,minimize=true)'
    catch
        println("percent must be between 0 and 1")
    end
    
    # Save posterior mean
    cov_filename = joinpath(outpath(m),"cov.h5")
    posterior_fid = h5open(cov_filename,"w") do posterior_fid
        posterior_fid["Θ_hat"] = convert(Matrix{Float32}, Θ_hat)
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
    
    # Open and start the file
    mainParams_out = joinpath(tablepath(m), "moments_mainParams.tex")
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
        if(ismatch(r"rho_chi",param.texLabel)) # ??? || (isequal(subspec,7) && texLabel == ":ρ_b"))
            continue
        end

        @printf(mainParams_fid, "\\\\ \n \$\%4.99s\$ & ", param.texLabel)
        
        #Print the values in outmat
        for val in outmat[index,:]
            @printf(mainParams_fid, "\%8.3f & ",val)
        end

        important_para = [important_para, index]
        
    end

    # Close the file
    endTexTableDoc(mainParams_fid;small=true)

    ###########################################################################################
    ## 4b. Write to Table 2: Prior mean, std dev and posterior mean, bands for other params
    ###########################################################################################
    
    periphParams_out = joinpath(tablepath(m), "moments_periphParams_0.tex")
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
            periphParams_out = joinpath(tablepath(m),filename)
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
    prioPostMean_out = joinpath(tablepath(m), "moments_prioPostMean.tex")
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
            prioPostMean_out = joinpath(tablepath(m),filename)
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
    println("Tables are in ",tablepath(m))

end


function find_density_bands(draws::Matrix, percent::Real; minimize::Bool=true)
    ## Returns a [2 x cols(draws)] matrix `bands` such that `percent` of the mass of `draws[:,i]` is above
    ## `bands[1,i]` and below `bands[2,i]`.
    ## 
    ## Inputs:
    ## -draws: [num_draws x num_draw_dimensions] matrix
    ## -percent: percent of data within bands (e.g. .9 to get 90% of mass within bands)
    ## -minimize: if =1, choose shortest interval, otherwise just chop off lowest and highest (percent/2)
    ##
    ## Output:
    ## -[2 x num_draw_dimensions] matrix 

    if(percent < 0 || percent > 1)
       throw(DomainError())
    end
    
    num_draws, num_draw_dimensions = size(draws)
    band  = zeros(2, num_draw_dimensions)
    num_in_band  = round(percent * num_draws)
    
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

function beginTexTableDoc(fid::IOStream)

    @printf(fid,"\\documentclass[12pt]{article}\n")
    @printf(fid,"\\usepackage[dvips]{color}\n")
    @printf(fid,"\\begin{document}\n")
    @printf(fid,"\\pagestyle{empty}\n")
    @printf(fid,"\\begin{table}[h] \\centering\n")
    
end

# Prints the necessarily lines to end a table and close a document
# to file descriptor fid and closes the file
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


function find_density_bands(draws::Matrix, percent::Real; minimize::Bool=true)
    ## Returns a [2 x cols(draws)] matrix `bands` such that `percent` of the mass of `draws[:,i]` is above
    ## `bands[1,i]` and below `bands[2,i]`.
    ## 
    ## Inputs:
    ## -draws: [num_draws x num_draw_dimensions] matrix
    ## -percent: percent of data within bands (e.g. .9 to get 90% of mass within bands)
    ## -minimize: if =1, choose shortest interval, otherwise just chop off lowest and highest (percent/2)
    ##
    ## Output:
    ## -[2 x num_draw_dimensions] matrix 

    if(percent < 0 || percent > 1)
       throw(DomainError())
    end
    
    num_draws, num_draw_dimensions = size(draws)
    band  = zeros(2, num_draw_dimensions)
    num_in_band  = round(percent * num_draws)
    
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

function beginTexTableDoc(fid::IOStream)

    @printf(fid,"\\documentclass[12pt]{article}\n")
    @printf(fid,"\\usepackage[dvips]{color}\n")
    @printf(fid,"\\begin{document}\n")
    @printf(fid,"\\pagestyle{empty}\n")
    @printf(fid,"\\begin{table}[h] \\centering\n")
    
end

# Prints the necessarily lines to end a table and close a document
# to file descriptor fid and closes the file
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
