"""
```
estimate_bma(m, data, prior; return_output = false, filestring_addl = [])
```

Estimate a Bayesian Model Average.

### Arguments:
- `m::PoolModel`: PoolModel object

### Optional Arguments:
- `data`: well-formed data as `Matrix` or `DataFrame`. If this is not provided, the `load_data` routine will be executed.
- `prior::Float64`: prior probability between the two models

### Keyword Arguments:
- `return_output::Bool`: option to return output. If false, `nothing` is returned.
- `filestring_addl::Vector{String}`: Additional strings to append to output files.
"""
function estimate_bma(m::PoolModel, df::DataFrame, prior::Float64 = 0.5;
                      return_output::Bool = false, save_output::Bool = true,
                      filestring_addl::Vector{String} = Vector{String}(undef,0))
    return estimate_bma(m, df_to_matrix(m, df), prior;
                        return_output = return_output, save_output = save_output,
                        filestring_addl = filestring_addl)
end

function estimate_bma(m::PoolModel, data::Matrix{Float64} = Matrix{Float64}(0,0),
                      prior::Float64 = 0.5; return_output::Bool = false, save_output::Bool = true,
                      filestring_addl::Vector{String} = Vector{String}(undef,0))

    # Compute λ weights
    λ = zeros(size(data,2))
    λ[1] = prior * data[1,1] / (prior * data[1,1] + (1 - prior) * data[2,1])
    for t in 2:size(data,2)
        λ[t] = λ[t-1] * data[1,t] / (λ[t-1] * data[1,t] + (1 - λ[t-1]) * data[2,t])
    end

    # Save weights
    if save_output
        h5open(rawpath(m, "estimate", "bmasave.h5", filestring_addl), "w") do file
            write(file, "bmaparams", λ)
        end
    end

    if return_output
        prob_vec = λ .* vec(data[1,:]) + (1 .- λ) .* vec(data[2,:])
        return λ, prob_vec
    end
end
