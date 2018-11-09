# Think about updating the reference trial writing
# to write reference trials to directories with the name of the git commit
# that the trial corresponds to.
# I.e. if the trial was run under commit f85e8e7, then within the top-level
# "ref" directory, there will be a directories named after commits
# Will also need to think about best ways to parse these directories
# in some sort of chronological order, or having a bunch of commit id labeled
# directories could get a bit unwieldy

using JLD2

"""
```
write_ref_trial(trial, trial_name)
```

Write a reference trial to a JLD file, to act as the standard that new trials are
benchmarked against.

### Arguments
- `trial::BenchmarkTools.Trial`: The trial object that is being written.
- `trial_name::String`: The name of the trial being written.

"""
function write_ref_trial(trial::BenchmarkTools.Trial, trial_name::String)
    filepath = "../reference/$trial_name.jld"
    if isfile(filepath)
        println("There is already a reference file in $filepath. Do you want to overwrite it? (y/n)")
        response = readline(STDIN)
        if response == "y"
            rm(filepath)
            op = "Overwrote"
        elseif response == "n"
            println("File will be left as is. Aborting write.")
            return nothing
        else
            throw("Invalid response, must be y or n")
        end
    else
        op = "Wrote"
    end

    d = Dict{Symbol, Union{Vector, Int}}()
    d[:times] = trial.times
    d[:gctimes] = trial.gctimes
    d[:memory] = trial.memory
    d[:allocs] = trial.allocs

    save(filepath, trial_name, d)

    println("$op "*filepath)
end

function write_ref_trial_group(group::Dict{Symbol, Dict},
                               trial_names::Vector{Symbol},
                               group_name::String)
    filepath = "../reference/$group_name.jld"
    if isfile(filepath)
        println("There is already a reference file in $filepath. Do you want to overwrite it? (y/n)")
        response = readline(STDIN)
        if response == "y"
            rm(filepath)
            op = "Overwrote"
        elseif response == "n"
            println("File will be left as is. Aborting write.")
            return nothing
        else
            throw("Invalid response, must be y or n")
        end
    else
        op = "Wrote"
    end

    save(filepath, group_name, group)
    println("$op "*filepath)
end
