using Documenter, DSGE, Distributions, DataFrames

makedocs(modules = [DSGE],
         clean = false,
         format = :html,
         sitename = "DSGE.jl",
         authors = "FRBNY-DSGE",
         linkcheck = !("skiplinks" in ARGS),
         strict = true
         pages = Any[
             "Home" => "index.md",
         ],
         doctest=  false # for now
)


deploydocs(
     repo = "github.com/FRBNY-DSGE/DSGE.jl",
     target = "build",
     deps = nothing
     make = nothing
)
