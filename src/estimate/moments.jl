##                           OVERVIEW
## moments.jl: Computes moments, tabulates parameter moments, and plots parameter
## draws from the prior and posterior distribution.

## IMPORTANT VARIABLES
## PLOTDRAWS: Specify whether or not you want to plot the parameter draws
## (runs momPlot.m)
## percent: In the output laTex table, we report <percent> bands for
##     parameter draws from the posterior

##     INPUTS
##     infile1: This refers to a an output file (mhparam*) from gibb.m,
##     which holds draws from the posterior.
##
##     infile4: This refers to a an output file (post*) from gibb.m, which
##     holds the value of the posterior, for each posterior draw.
##
##     theta: (ndraws x npara) matrix holding the posterior draws
##	   (stored in mhpara*, output from gibb.m)
##
##     post: (ndraws x 1) matrix holding the posterior value
##	   (stored in post*, output from gibb.m)

##     OUTPUTS
##     From makeMomentsTables():
##
##     1: (mhpara*Mean): holds the mean parameter across
##     draws from the posterior.
##
##     2: (*Mom_MainParams*): laTex table that lists the
##	   moments for important parameters.
##
##     3: (*Mom_PeriphParams*): laTex table that lists the
##	       moments for less important parameters.
##		   4: (*PrioPostMean*) that lists the prior and
##		   posterior means

##		  



using HDF5
using Debug


@debug function computeMoments{T<:AbstractDSGEModel}(m::T, percent::Float64 = 0.90)
    
    ## Computes prior and posterior parameter moments, tabulates them in various TeX tables,
    ## and plots parameter draws from the prior and posterior distribution
    ##
    ## Inputs:
    ## - m: an instance of an AbstractDSGEModel subtype
    ## - percent: the percentage of the mass of draws from Metropolis-Hastings
    ##            we want to include between bands shown in the table

    
    # Read in the matrix of parameter draws from metropolis-hastings

    infile = joinpath(outpath(),"sim_save.h5")
    println("Reading draws from Metropolis-Hastings from $infile...")

    Θ = []
    post = []

    try
        fid = h5open(infile,"r") do fid
            Θ = read(fid,"parasim")
            post = read(fid,"postsim")
        end
    catch
        @printf(1,"Could not open file %s", infile)
    end

    # Convert back to Float64 for compatability with other variables
    Θ = convert(Matrix{Float64},Θ)

    num_draws = size(Θ,1)

    # Produce TeX table of moments
    makeMomentTables(m,Θ,percent)

    # ?? Do we want to calc covariance here? Haven't we done it already in estimate?

    # ?? What to return?
    #return Θ, post
end




@debug function makeMomentTables{T<:FloatingPoint}(m::AbstractDSGEModel, Θ::Array{T,2}, percent::Float64)
    
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


    ## Step 1: Extract moments of prior distribution from m.parameters

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

            prior_means[i] = param.priordist.α
            prior_stddev[i] = param.priordist.β

        elseif isa(param.priordist, Distributions.Gamma)

            prior_means[i] = param.priordist.α
            prior_stddev[i] = param.priordist.θ  # small \theta
            
        end
    end
    
    
    ## Step 2: Compute moments and `percent' bands from parameter draws

    Θ_hat = mean(Θ,1)'                    # posterior mean for each parameter
    Θ_sig = cov(Θ, mean=Θ_hat')           # posterior std dev for each parameter
    Θ_bands = []
    
    # ?? Do we want to be doing error handling like this?
    try
        Θ_bands = find_density_bands(Θ,percent,minimize=true)'
    catch
        println("percent must be between 0 and 1")
    end
    
    # Save posterior mean
    cov_filename = joinpath(outpath(),"cov.h5")
    posterior_fid = h5open(cov_filename,"w") do posterior_fid
        posterior_fid["Θ_hat"] = convert(Matrix{Float32}, Θ_hat)
    end


    
    ## Step 3: Create variables for writing to output table
    colnames = ["Parameter ", "Prior Mean ", "Prior Stdd ", "Post Mean ", "$(100*percent)\\% {\\tiny Lower Band} ", "$(100*percent)\\% {\\tiny Upper Band} " ]

    outmat = [prior_means prior_stddev Θ_hat Θ_bands]      # prior mean and std dev, posterior mean, and bands (n_params x 5)
    outmat2 = [prior_means Θ_hat]                   # prior mean and posterior mean (n_params x 2)

    
    ## Step 4: Write to Table 1: prior mean, std dev, posterior mean and bands for IMPORTANT parameters

    mainParams_out = joinpath(tablepath(), "moments_mainParams.tex")
    mainParams_fid = open(mainParams_out,"w")

    # Preamble
    @printf(mainParams_fid,"\\documentclass[12pt]{article}\n")
    @printf(mainParams_fid,"\\usepackage[dvips]{color}\n")
    @printf(mainParams_fid,"\\begin{document}\n")
    @printf(mainParams_fid,"\\pagestyle{empty}\n")
    @printf(mainParams_fid,"\\begin{table}[h] \\centering\n")
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
            !ismatch(r"psi", param.texLabel) &&
            !ismatch(r"nu_l", param.texLabel) &&
            !ismatch(r"pistar", param.texLabel) &&
            (!ismatch(r"ups", param.texLabel)))  ##Is this correct? mspec == 16 or u_^*
            continue
        end
            
        # TODO: Decide whether subspec should be a field in the model
        if(ismatch(r"\\rho_chi",param.texLabel)) # ??? || (isequal(subspec,7) && texLabel == ":ρ_b"))
            continue
        end

        @printf(mainParams_fid, "\\\\ \n \$\%4.99s\$ & ", param.texLabel)
        
        #Print the values in outmat
        for val in outmat[index,:]
            @printf(mainParams_fid, "\%8.3f & ",val)
        end

        important_para = [important_para, index]
        
    end

    @printf(mainParams_fid,"\\\\  \\\hline\n")
    @printf(mainParams_fid,"\\end{tabular}}\n")
    @printf(mainParams_fid,"\\end{table}\n")
    @printf(mainParams_fid,"\\end{document}")
    close(mainParams_fid)    
    
    # Write to Table 2: Prior mean, std dev and posterior mean, bands for other params
    periphParams_out = joinpath(tablepath(), "moments_periphParams_0.tex")
    periphParams_fid = open(periphParams_out,"w")
    
    @printf(periphParams_fid,"\\documentclass[12pt]{article}\n")
    @printf(periphParams_fid,"\\usepackage[dvips]{color}\n")
    @printf(periphParams_fid,"\\begin{document}\n")
    @printf(periphParams_fid,"\\pagestyle{empty}\n")
    @printf(periphParams_fid,"\\begin{table}[h] \\centering\n")
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
            
            @printf(periphParams_fid,"\\\\  \\\hline\n")
            @printf(periphParams_fid,"\\end{tabular}}\n")
            @printf(periphParams_fid,"\\end{table}\n")
            @printf(periphParams_fid,"\\end{document}")
            
            table_count = table_count + 1
            close(periphParams_fid)
            
            filename = @sprintf("moments_periphParams_%d.tex",table_count)
            periphParams_out = joinpath(tablepath(),filename)
            periphParams_fid = open(periphParams_out,"w")
            
            @printf(periphParams_fid,"\\documentclass[12pt]{article}\n")
            @printf(periphParams_fid,"\\usepackage[dvips]{color}\n")
            @printf(periphParams_fid,"\\begin{document}\n")
            @printf(periphParams_fid,"\\pagestyle{empty}\n")
            @printf(periphParams_fid,"\\begin{table}[h] \\centering\n")
            @printf(periphParams_fid,"\\caption{Parameter Estimates}\n")
            @printf(periphParams_fid,"\\vspace*{.2cm}\n")
            @printf(periphParams_fid,"{\\small \n")
            @printf(periphParams_fid,"\\begin{tabular}{lllllll}\\hline \n")
            
            for col in colnames
                @printf(periphParams_fid, "%4.99s & ", col)
            end
        end
        
        #print the parameter name
        @printf(periphParams_fid, "\\\\ \n \$\%4.99s\$ & ", param.texLabel)
        
        #Print the values in outmat
        for val in outmat[index,:]
            @printf(periphParams_fid, "\%8.3f & ",val)
        end
        
        other_para = other_para+1
        
   end
   
   @printf(periphParams_fid,"\\\\  \\\hline\n")
   @printf(periphParams_fid,"\\end{tabular}}\n")
   @printf(periphParams_fid,"\\end{table}\n")
   @printf(periphParams_fid,"\\end{document}")
   close(periphParams_fid)    


   ## Write to Table 3: Prior mean and posterior mean for all parameters

    prioPostMean_out = joinpath(tablepath(), "moments_prioPostMean.tex")
    prioPostMean_fid = open(prioPostMean_out,"w")


    # Counter for parameters to track length of table and number of tables in excess of
    # default (1)
    table_count = 0

    @printf(prioPostMean_fid,"\\documentclass[12pt]{article}\n")
    @printf(prioPostMean_fid,"\\usepackage[dvips]{color}\n")
    @printf(prioPostMean_fid,"\\begin{document}\n")
    @printf(prioPostMean_fid,"\\pagestyle{empty}\n")
    @printf(prioPostMean_fid,"\\begin{table}[h] \\centering\n")
    @printf(prioPostMean_fid,"\\caption{Parameter Estimates: Prior and Posterior Mean}\n")
    @printf(prioPostMean_fid,"\\vspace*{.5cm}\n")
    @printf(prioPostMean_fid,"\\begin{tabular}{ccc}\\hline \n")
    @printf(prioPostMean_fid," Parameter & Prior & Posterior  \\tabularnewline \\hline\n")
    
    for (index, param) in enumerate(m.parameters)

        if ((index % 40 == 0) && (index != length(m.parameters)) )

            @printf(prioPostMean_fid,"\\\\  \\\hline\n")
            @printf(prioPostMean_fid,"\\end{tabular}}\n")
            @printf(prioPostMean_fid,"\\end{table}\n")
            @printf(prioPostMean_fid,"\\end{document}")
            
            table_count = table_count + 1
            close(prioPostMean_fid)
            
            filename = @sprintf("moments_prioPostMean_%d.tex",table_count)
            prioPostMean_out = joinpath(tablepath(),filename)
            prioPostMean_fid = open(prioPostMean_out,"w")
            
            @printf(prioPostMean_fid,"\\documentclass[12pt]{article}\n")
            @printf(prioPostMean_fid,"\\usepackage[dvips]{color}\n")
            @printf(prioPostMean_fid,"\\begin{document}\n")
            @printf(prioPostMean_fid,"\\pagestyle{empty}\n")
            @printf(prioPostMean_fid,"\\begin{table}[h] \\centering\n")
            @printf(prioPostMean_fid,"\\caption{Parameter Estimates: Prior and Posterior Mean}\n")
            @printf(prioPostMean_fid,"\\vspace*{.5cm}\n")
            @printf(prioPostMean_fid,"\\begin{tabular}{ccc}\\hline \n")
            @printf(prioPostMean_fid," Parameter & Prior & Posterior  \\tabularnewline \\hline\n")

        end

        
        @printf(prioPostMean_fid, "\\\\ \n \$\%4.99s\$ & ", param.texLabel)

        val = outmat2[index,:]
        @printf(prioPostMean_fid, "\%8.3f &   \%8.3f  ", val[1], val[2])
        
    end

    @printf(prioPostMean_fid, "\\\\ \\\hline\n")
    @printf(prioPostMean_fid,"\\end{tabular}\n")
    @printf(prioPostMean_fid,"\\end{table}\n")
    @printf(prioPostMean_fid,"\\end{document}")
    close(prioPostMean_fid)

   println("Tables are in ",tablepath())

end


## function plotParamDraws()
    
## end

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

