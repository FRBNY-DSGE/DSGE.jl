using Pkg
Pkg.add("Documenter")
using Documenter, DSGE

makedocs(modules = [DSGE],
         clean = false,
         format = Documenter.HTML(),
         sitename = "DSGE.jl",
         authors = "FRBNY-DSGE",
         linkcheck = false,
         strict = false,
         pages = Any[
                     "Home"                                   => "index.md",
                     "Learning How to Use DSGE.jl"            => "learning_how_to_use_dsgejl.md",
                     "Model Design"                           => "model_design.md",
                     "Special Model Types"                    => "special_model_types.md",
                     "Model Implementation Details"           => "model_implementation_details.md",
                     "Running An Existing Model"              => "running_existing_model.md",
                     "Input Data"                             => "input_data.md",
                     "FRBNY Model Input Data"                 => "frbny_data.md",
                     "Solving the Model"                      => "solving.md",
                     "Estimation"                             => "estimation.md",
                     "Forecasting"                            => "forecast.md",
                     "Impulse Response Functions"             => "irf.md",
                     "Computing Means and Bands"              => "means_bands.md",
                     "Alternative Policies"                   => "altpolicy.md",
                     "Alternative Scenarios"                  => "scenarios.md",
                     "Forecast Decomposition"                 => "forecast_decomposition.md",
                     "Plotting"                               => "plotting.md",
                     "Algorithms"                             => "algorithms.md",
                     "Advanced Usage"                         => "advanced_usage.md",
                     "Contributing to DSGE.jl"                => "contributing.md",
                     "MATLAB to Julia Transition: Estimation" => "MatlabToJuliaTransition.md",
                     "MATLAB to Julia Transition: Forecast"   => "julia_forecasting.md",
                     "License"                                => "license.md"
         ],
         doctest = false # for now
)

# if "deploy" in ARGS
#    include("../../fake_travis.jl")
# end
deploydocs(
    repo = "github.com/FRBNY-DSGE/DSGE.jl.git",
    target = "build",
    deps = nothing,
    devbranch = "main",
    branch = "gh-pages",
    # versions = "v#",
    # julia = "0.7",
    # osname = "osx",
    make = nothing
)
