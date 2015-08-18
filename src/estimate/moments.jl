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
##	   priotheta: (ndraws x npara) matrix holding the prior draws
##	   (stored in mhNIpriopara*, output from prio.m)
##	   (optional input)

##     OUTPUTS
##     From momTable.m:
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

## TODO: initializePrograms etc

# Computes moments and outputs them as a TeX table
function computeMoments(verbose::Bool = 0, percent::Float64 = 0.90)

    ## Input file

    infile = joinpath("$savepath","sim_save.h5")
    println("Infile: $infile")

    θ = []
    θ_prior = []
    post = []

    try
        fid = h5open(infile,"r") do fid
            θ = read(fid,"parasim")
            post = read(fid,"postsim")
        end
    catch
        @printf "Could not open file #infile"
    end

    nDraws = size(θ,1)
    
    ## Produce TeX table of moments
    makeMomentsTable(θ)
    
end

## ?? Float32s
function makeMomentsTable(θ::Vector{Param})
    
    # Variables we have to somehow set:
    percent = 0.4  # I made this up
    pmean = rand(spec["n_params"])
    
    # ?? Do we have to do this for the prior or is it done in something like initializePrograms?
    # Compute mean, std dev, and 90% bands across draws from the posterior
    
    θ_hat = mean(θ,1)'       # posterior mean
    θ_sig = cov(θ,1)         # posterior std dev

    # ??? what's hpdint in MATLAB doing?
    θ_bands = []

    # Save posterior mean
    posterior_fid = h5open("$outpath","w") do posterior_fid
        posterior_fid["θ_hat"] = convert(Matrix{Float32}, θ_hat)
    end

    ## Open tex files
    mainParams_out = joinpath("$tablepath", "moments_mainParams.tex")
    mainParams_fid = open(mainParams_out,"w")

    periphParams_out = joinpath("$tablepath", "moments_periphParams.tex")
    periphParams_fid = open(periphParams_out,"w")

    prioPostMean_out = joinpath("$tablepath", "moments_prioPostMean.tex")
    prioPostMean_fid = open(prioPostMean_out,"w")

    ## Create variables for writing to output table
    colnames = ["Parameter ", "Prior Mean ", "Prior Stdd ", "Post Mean ", "$(100*percent)\\% {\\tiny Lower Band} ", "$(100*percent)\\% {\\tiny Upper Band} " ]

    outmat = [pmean, pstdd, θ_hat, θ_bands]       # prior mean and std dev, posterior mean, and bands
    outmat2 = [pmean, θ_hat]                   # prior mean and posterior mean

    ## Write to Table 1: prior mean, std dev, posterior mean and bands for IMPORTANT parameters
    
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
    
    
    ## Keep track of indices for important parameters within para_names
    important_para = [];
    
    # Go through all params in model, check tex names
    ## for param in model.Θ
        
    ##     texLabel = param.texLabel
        
    ##     # ?? rho_ might need to be ρ
    ##     if (!ismatch(r"rho_", texLabel) &&
    ##         !ismatch(r"zeta_", texLabel) &&
    ##         !ismatch(r"psi_", texLabel) &&
    ##         !ismatch(r"nu_l", texLabel) &&
    ##         !ismatch(r"pi^\*_", texLabel) &&
    ##         ( model.spec != 16 || !ismatch(r"u^\*", texLabel))) 
    ##         continue
    ##     end
            
    ##     # !! What is subspec in the Julia code?
    ##     if(texLabel == "\rho_{\chi}" || (isequal(subspec,7) && texLabel == "\rho_{b}"))
    ##         continue
    ##     end

    ##     # ?? figure out how to print a whole row of a matrix! 
    ##     for i in outmat
    ##         @printf(mainParams_fid,"\\\\ \\n %4.99s &  ", texLabel)
    ##         @printf(mainParams_fid, "8.3f & ", outmat[param,i])
    ##     end

    ##     important_para = [important_para, param]
        
    ## end
    
    @printf(mainParams_fid,"\\\\  \\\hline\n")
    @printf(mainParams_fid,"\\end{tabular}}\n")
    @printf(mainParams_fid,"\\end{table}\n")
    @printf(mainParams_fid,"\\end{document}")
    close(mainParams_fid)    
    
    # Write to Table 2: Prior mean, std dev and posterior mean, bands for other params
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
        @printf(periphParams_fid, "\$\%4.99s\$ & ", col)
    end
    
    # Counter for parameters to track length of table and number of tables in excess of
    # default (1)
    other_para = 1
    table_count = 0
    for param in Θ
        if in(param, important_para)
            continue
        end
            
        if (other_para%25 == 0) && !(i==length(Θ))
            @printf(periphParams_fid,"\\\\  \\\hline\n")
            @printf(periphParams_fid,"\\end{tabular}}\n")
        @printf(periphParams_fid,"\\end{table}\n")
        @printf(periphParams_fid,"\\end{document}")
        table_count = table_count + 1
        close(periphParams_fid)
        filename = @sprintf("periphParams_%d.tex",table_count)
        periphParams_out = joinpath("$tablepath",filename)
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
                @printf(periphParams_fid, "\$%4.99s\$ & ", col)
            end
        
        
        end
        @printf(periphParams_fid, "\\\\ \n \$\%4.99s\$ & ",param)
        for val in outmat[param,:]
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
   
    # Column names
    for col in colnames
        @printf(mainParams_fid, "%4.99s & ", col)
    end
    
    
    ## Keep track of indices for important parameters within para_names
    important_para = [];
    
    # Go through all params in model, check tex names
    ## for param in model.Θ
        
    ##     texLabel = param.texLabel
        
    ##     # ?? rho_ might need to be ρ
    ##     if (!ismatch(r"rho_", texLabel) &&
    ##         !ismatch(r"zeta_", texLabel) &&
    ##         !ismatch(r"psi_", texLabel) &&
    ##         !ismatch(r"nu_l", texLabel) &&
    ##         !ismatch(r"pi^\*_", texLabel) &&
    ##         ( model.spec != 16 || !ismatch(r"u^\*", texLabel))) 
    ##         continue
    ##     end
            
    ##     # !! What is subspec in the Julia code?
    ##     if(texLabel == "\rho_{\chi}" || (isequal(subspec,7) && texLabel == "\rho_{b}"))
    ##         continue
    ##     end

    ##     # ?? figure out how to print a whole row of a matrix! 
    ##     for i in outmat
    ##         @printf(mainParams_fid,"\\\\ \\n %4.99s &  ", texLabel)
    ##         @printf(mainParams_fid, "8.3f & ", outmat[param,i])
    ##     end

    ##     important_para = [important_para, param]
        
    ## end
    
    @printf(mainParams_fid,"\\\\  \\\hline\n")
    @printf(mainParams_fid,"\\end{tabular}}\n")
    @printf(mainParams_fid,"\\end{table}\n")
    @printf(mainParams_fid,"\\end{document}")
    close(mainParams_fid)    
    
    # Write to Table 2: Prior mean, std dev and posterior mean, bands for other params
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
        @printf(periphParams_fid, "\$\%4.99s\$ & ", col)
    end
    
    # Counter for parameters to track length of table and number of tables in excess of
    # default (1)
    other_para = 1
    table_count = 0
    for param in Θ
        if in(param, important_para)
            continue
        end
            
        if (other_para%25 == 0) && !(i==length(Θ))
            @printf(periphParams_fid,"\\\\  \\\hline\n")
            @printf(periphParams_fid,"\\end{tabular}}\n")
        @printf(periphParams_fid,"\\end{table}\n")
        @printf(periphParams_fid,"\\end{document}")
        table_count = table_count + 1
        close(periphParams_fid)
        filename = @sprintf("periphParams_%d.tex",table_count)
        periphParams_out = joinpath("$tablepath",filename)
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
                @printf(periphParams_fid, "\$%4.99s\$ & ", col)
            end
        
        
        end
        @printf(periphParams_fid, "\\\\ \n \$\%4.99s\$ & ",param)
        for val in outmat[param,:]
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


    # Column names
    for col in colnames
        @printf(mainParams_fid, "%4.99s & ", col)
    end
    
    
    ## Keep track of indices for important parameters within para_names
    important_para = [];
    
    # Go through all params in model, check tex names
    ## for param in model.Θ
        
    ##     texLabel = param.texLabel
        
    ##     # ?? rho_ might need to be ρ
    ##     if (!ismatch(r"rho_", texLabel) &&
    ##         !ismatch(r"zeta_", texLabel) &&
    ##         !ismatch(r"psi_", texLabel) &&
    ##         !ismatch(r"nu_l", texLabel) &&
    ##         !ismatch(r"pi^\*_", texLabel) &&
    ##         ( model.spec != 16 || !ismatch(r"u^\*", texLabel))) 
    ##         continue
    ##     end
            
    ##     # !! What is subspec in the Julia code?
    ##     if(texLabel == "\rho_{\chi}" || (isequal(subspec,7) && texLabel == "\rho_{b}"))
    ##         continue
    ##     end

    ##     # ?? figure out how to print a whole row of a matrix! 
    ##     for i in outmat
    ##         @printf(mainParams_fid,"\\\\ \\n %4.99s &  ", texLabel)
    ##         @printf(mainParams_fid, "8.3f & ", outmat[param,i])
    ##     end

    ##     important_para = [important_para, param]
        
    ## end
    
    @printf(mainParams_fid,"\\\\  \\\hline\n")
    @printf(mainParams_fid,"\\end{tabular}}\n")
    @printf(mainParams_fid,"\\end{table}\n")
    @printf(mainParams_fid,"\\end{document}")
    close(mainParams_fid)    
    
    # Write to Table 2: Prior mean, std dev and posterior mean, bands for other params
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
        @printf(periphParams_fid, "\$\%4.99s\$ & ", col)
    end
    
    # Counter for parameters to track length of table and number of tables in excess of
    # default (1)
    other_para = 1
    table_count = 0
    for param in Θ
        if in(param, important_para)
            continue
        end
            
        if (other_para%25 == 0) && !(i==length(Θ))
            @printf(periphParams_fid,"\\\\  \\\hline\n")
            @printf(periphParams_fid,"\\end{tabular}}\n")
        @printf(periphParams_fid,"\\end{table}\n")
        @printf(periphParams_fid,"\\end{document}")
        table_count = table_count + 1
        close(periphParams_fid)
        filename = @sprintf("periphParams_%d.tex",table_count)
        periphParams_out = joinpath("$tablepath",filename)
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
                @printf(periphParams_fid, "\$%4.99s\$ & ", col)
            end
        
        
        end
        @printf(periphParams_fid, "\\\\ \n \$\%4.99s\$ & ",param)
        for val in outmat[param,:]
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

    @printf(prioPostMean_fid,"\\documentclass[12pt]{article}\n")
    @printf(prioPostMean_fid,"\\usepackage[dvips]{color}\n")
    @printf(prioPostMean_fid,"\\begin{document}\n")
    @printf(prioPostMean_fid,"\\pagestyle{empty}\n")
    @printf(prioPostMean_fid,"\\begin{table}[h] \\centering\n")
    @printf(prioPostMean_fid,"\\caption{Parameter Estimates: Prior and Posterior Mean}\n")
    @printf(prioPostMean_fid,"\\vspace*{.5cm}\n")
    @printf(prioPostMean_fid,"\\begin{tabular}{ccc}\\hline \n")
    @printf(prioPostMean_fid," Parameter & Prior \n")
    @printf(prioPostMean_fid,"\\\\  \\\hline\n")
    for param in Θ
        @printf(prioPostMean_fid, " \\n \$%4.99s\$ ", param )
        for val in outmat2[param,:]
            @printf(periphParams_fid, "\%8.3f & ",val)
        end
        @printf(prioPostMean_fid, "\\\\ ")
    end
    @printf(prioPostMean_fid, "\\\\ \\hline")
    close(prioPostMean_fid)
    @printf("Tables are in $tablepath")
end


function plotParamDraws()

end
