using BenchmarkTools

function calculate_color(diff::Float64, title::Symbol)
    threshold = improvement_thresholds(title)
    if diff > threshold
        color = :red
    elseif diff < -threshold
        color = :green
    else
        color = :normal
    end
    return color
end

function trial_to_dict(trial::BenchmarkTools.Trial)
    d = Dict{Symbol, Union{Vector, Int}}()
    d[:times] = trial.times
    d[:gctimes] = trial.gctimes
    d[:memory] = trial.memory
    d[:allocs] = trial.allocs

    return d
end

function construct_trial_group(trials::Vector{BenchmarkTools.Trial},
                               component_names::Vector{Symbol})
    d = Dict{Symbol, Dict}()
    for (trial, component) in zip(trials, component_names)
        d[component] = trial_to_dict(trial)
    end
    return d
end

function improvement_thresholds(title::Symbol)
    title_str = string(title)
    if contains(title_str, "gctimes")
        threshold = 100.
    elseif contains(title_str, "memory")
        threshold = 5.
    elseif contains(title_str, "allocs")
        threshold = 5.
    else # if the title is for times
        threshold = 10.
    end
    return threshold
end

function stat_maps()
    d = Dict()
    d[:min] = minimum
    d[:median] = median
    d[:mean] = mean
    d[:max] = maximum
    d[:identity] = identity

    return d
end
