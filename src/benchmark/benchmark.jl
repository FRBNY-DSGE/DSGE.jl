using BenchmarkTools
using JLD2

"""
```
print_all_benchmarks(trial, ref, trial_name; stats = [:min, :median, :mean, :max], factors = [:times, :gctimes, :memory, :allocs])
```

For printing the benchmark comparisons given by the BenchmarkTools package of a single trial.

### Arguments
- `trial::BenchmarkTools.Trial`: A trial of type BenchmarkTools.Trial, outputted from `@benchmark`.
- `ref::String`: The location of the previous trial, saved by the `write_ref_trial` method
- `trial_name::String`: The name of the trial being benchmarked, i.e. `forecast_one`, if the `forecast_one` function is the trial being benchmarked
"""
function print_all_benchmarks(trial::BenchmarkTools.Trial, ref::String,
                              trial_name::String;
                              stats::Vector{Symbol} = [:min, :median, :mean, :max],
                              factors::Vector{Symbol} = [:times, :gctimes, :memory, :allocs])
    print_all_benchmarks(trial_to_dict(trial), ref, trial_name, stats = stats, factors = factors)
end

# For a single trial comparison
function print_all_benchmarks(trial::Dict{Symbol, Union{Vector, Int}}, ref::String,
                              trial_name::String;
                              stats::Vector{Symbol} = [:min, :median, :mean, :max],
                              factors::Vector{Symbol} = [:times, :gctimes, :memory, :allocs])
    prev_trial = load(ref, trial_name)
    print_all_benchmarks(trial, prev_trial, trial_name, stats = stats, factors = factors)
end


"""
```
print_all_benchmarks(trial, prev_trial, trial_name; stats = [:min, :median, :mean, :max], factors = [:times, :gctimes, :memory, :allocs])
```

The underlying method for printing the benchmark comparisons given by the BenchmarkTools package of a single trial.

### Arguments
- `trial::Dict{Symbol, Union{Vector, Int}}`: A dictionary mapping the `factors` to their values for an individual trial.
- `prev_trial::Dict{Symbol, Union{Vector, Int}}`: The comparison factors that `trial` is being benchmarked against.
- `trial_name::String`: The name of the trial that is being benchmarked, i.e. `forecast_one`, if the `forecast_one` function is the trial being benchmarked.

### Keyword Arguments
- `stats::Vector{Symbol}`: The statistics that can be shown when printing the benchmarks, which include `[:min, :median, :mean, :max]`.
- `factors::Vector{Symbol}`: The different factors within the trial that can be shown when printing the benchmarks, which include `[:times, :gctimes, :memory, :allocs]`.
- `group::Bool`: Whether or not this benchmark was part of a group of trials being benchmarked, for the purposes of indentation of the stdout.
"""
function print_all_benchmarks(trial::Dict{Symbol, Union{Vector, Int}},
                              prev_trial::Dict{Symbol, Union{Vector, Int}},
                              trial_name::String;
                              stats::Vector{Symbol} = [:min, :median, :mean, :max],
                              factors::Vector{Symbol} = [:times, :gctimes, :memory, :allocs],
                              group::Bool = false)
    # Print formatting of the trial sections
    trial_delimiter = repeat("-", length(trial_name))
    trial_name = group ? "  "*trial_name : trial_name
    trial_delimiter = group ? "  "*trial_delimiter : trial_delimiter
    println(trial_delimiter*"\n")
    print_with_color(:bold, trial_name)
    print("\n")

    # Printing benchmarks for each factor within a single trial
    for f in factors
        times = (f == :times)
        if f in [:memory, :allocs]
            print_single_benchmark(trial, prev_trial, :identity, f, times, group)
            print("\n")
            continue
        end
        for s in stats
            print_single_benchmark(trial, prev_trial, s, f, times, group)
        end
        print("\n")
    end
end

"""
```
print_all_benchmarks(group, ref, group_name, trial_names; stats = [:min, :median, :mean, :max], factors = [:times, :gctimes, :memory, :allocs])
```

For printing the benchmark comparisons given by the BenchmarkTools package of a group of trials.

### Arguments
- `group::Dict{Symbol, Dict}`: A dictionary mapping the `trial_names` to their corresponding trials, which are of type, `Dict{Symbol, Union{Vector, Int}}`.
- `ref::String`: The location of the .jld file containing the previous group to be benchmarked against.
- `trial_names::Vector{Symbol}`: The names of the trials being benchmarked, i.e. `forecast_one`, if the `forecast_one` function is one of the trials being benchmarked.

### Keyword Arguments
- `stats::Vector{Symbol}`: The statistics that can be shown when printing the benchmarks, which include `[:min, :median, :mean, :max]`.
- `factors::Vector{Symbol}`: The different factors within the trial that can be shown when printing the benchmarks, which include `[:times, :gctimes, :memory, :allocs]`.
"""
function print_all_benchmarks(group::Dict{Symbol, Dict}, ref::String,
                              group_name::String,
                              trial_names::Vector{Symbol};
                              stats::Vector{Symbol} = [:min, :median, :mean, :max],
                              factors::Vector{Symbol} = [:times, :gctimes, :memory, :allocs])
    # Load in the previous group to benchmark against
    prev_group = load(ref, group_name)

    # Print the group_name section header
    print_with_color(:green, group_name, bold = true)
    print("\n")

    # For each trial within a group, print all of the benchmarks.
    for trial_name in trial_names
        print_all_benchmarks(group[trial_name], prev_group[trial_name], string(trial_name),
                             stats = stats, factors = factors, group = true)
    end
end

"""
```
print_single_benchmark(trial, prev_trial, stat, factor, times)
```

### Arguments
- `trial::Dict{Symbol, Union{Vector, Int}}`: A dictionary mapping the `factors` to their values for an individual trial.
- `prev_trial::Dict{Symbol, Union{Vector, Int}}`: The comparison factors that `trial` is being benchmarked against.
- `stats::Symbol`: The statistic that will be shown when printing the benchmarks, one of `[:min, :median, :mean, :max]`.
- `factors::Symbol`: The factor within the trial that will be shown when printing the benchmarks, one of `[:times, :gctimes, :memory, :allocs]`.
- `times::Bool`: Whether or not the factor chosen was the `:times` factor, for the purposes of printing more grammatically correct output.
- `group::Bool`: Whether or not this benchmark was part of a group of trials being benchmarked, for the purposes of indentation of the stdout.
"""
function print_single_benchmark(trial::Dict{Symbol, Union{Vector, Int}},
                                prev_trial::Dict{Symbol, Union{Vector, Int}},
                                stat::Symbol, factor::Symbol,
                                times::Bool = false,
                                group::Bool = false)
    title = stat == :identity ? factor : Symbol(stat, " ", factor)
    calc_stat(x) = stat_maps()[stat](x)
    print_single_benchmark(title, calc_stat(trial[factor]),
                           calc_stat(prev_trial[factor]), times, group)
end

function print_single_benchmark(title::Symbol, new::T, old::T, times::Bool = false,
                                group::Bool = false) where T<:Real
    bad_state = times ? "slower" : "larger"
    good_state = times ? "faster" : "smaller"

    diff = round((new/old - 1)*100, 2)
    if diff > 0
        state = bad_state
    else
        state = good_state
    end
    color = calculate_color(diff, title)

    # Make sure the printing is indented properly to more aesthetically display group trial results
    if group
        print_with_color(color, "  $title is $(abs(diff))% $state\n")
    else
        print_with_color(color, "$title is $(abs(diff))% $state\n")
    end
end
