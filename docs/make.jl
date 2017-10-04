using Documenter, DSGE, Distributions, DataFrames

makedocs(modules = [DSGE],
         clean = false,
         format = :html,
         sitename = "DSGE.jl",
         authors = "FRBNY-DSGE",
         linkcheck = false,
         strict = false,
         pages = Any[
                     "Home"                                   => "index.md",
                     "Model Design"                           => "model_design.md",
                     "Running An Existing Model"              => "running_existing_model.md",
                     "Advanced Usage"                         => "advanced_usage.md",
                     "Input Data"                             => "input_data.md",
                     "FRBNY Model Input Data"                 => "frbny_data.md",
                     "Implementation Details"                 => "implementation_details.md",
                     "Solving the Model"                      => "solving.md",
                     "Estimating the Model"                   => "estimation.md",
                     "Forecasting"                            => "forecast.md",
                     "Computing Means and Bands"              => "means_bands.md",
                     "Alternative Policies"                   => "altpolicy.md",
                     "Alternative Scenarios"                  => "scenarios.md",
                     "Plotting"                               => "plotting.md",
                     "Algorithms"                             => "algorithms.md",
                     "Contributing to DSGE.jl"                => "contributing.md",
                     "MATLAB to Julia Transition: Estimation" => "MatlabToJuliaTransition.md",
                     "MATLAB to Julia Transition: Forecast"   => "julia_forecasting.md",
                     "License"                                => "license.md"
         ],
         doctest = false # for now
)


deploydocs(
     repo = "github.com/FRBNY-DSGE/DSGE.jl.git",
     target = "build",
     deps = nothing,
     julia = "0.6",
     osname = "osx",
     make = nothing
)