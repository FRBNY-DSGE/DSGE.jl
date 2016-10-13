using Documenter, DSGE, Distributions, DataFrames

makedocs(modules = [DSGE],
         clean = false,
         format = Documenter.Formats.HTML,
         sitename = "DSGE.jl",
         authors = "FRBNY-DSGE",
         linkcheck = !("skiplinks" in ARGS),
         strict = true,
         pages = Any[
                     "Home" => "intro.md",
                     "Model Design" => "model_design.md",
                     "Replicating the FRBNY DSGE Model" => "replicating_frbny.md",
                     "Advanced Usage" => "advanced_usage.md",
                     "Input Data" => "input_data.md",
                     "Implementation Details" => "implementation_details.md",
                     "Full API" => "index.md"
         ],
         doctest = false # for now
)


deploydocs(
     repo = "github.com/FRBNY-DSGE/DSGE.jl.git",
     target = "build",
     deps = nothing,
     julia = "0.4",
     osname = "osx",      
     make = nothing
)
