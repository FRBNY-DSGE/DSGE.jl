function load_cond_data(m, cond_type)
    # Preapare file name
    cond_vintage = get_setting(m, :cond_vintage)
    cond_id = get_setting(m, :cond_id)
    cond_file_name = inpath(m, "cond", "cond_vint=$(cond_vintage)_cdid=$(cond_id).csv")

    # Read data
    cond_data_raw, cond_headers = readcsv(cond_file_name; header=true)
    cond_data = nan(size(cond_data_raw))

    # Rearrange columns of conditional data to match indices defined in m.observables
    for (i,header) in enumerate(cond_headers)
        # If this column name is found in the m.observables dictionary, then move the column
        # to the correct index in the data matrix.
        key = symbol(header)
        if haskey(m.observables, key)
            j = m.observables[key]
            cond_data[:,j] = cond_data_raw[:,i]
        else
            # If this column is not found, determine nearest names for column and raise
            # warning.
            nearest_str = string_nearest(m.observables, header)
            warn("Invalid conditional data column: " * header) 
            warn("\tCandidates include: " * join(nearest_str, ", ") * ".")
        end
    end

    # If we are using semi-conditional data only, then we must index out the
    # semi-conditional columns and replace all other columns with NaN.
    if cond_type == :semi
        cond_semi_names = get_setting(m, :cond_semi_names)
        cond_semi_inds = [m.observables[sym] for sym in cond_semi_names]
        cond_semi_inds = setdiff(1:size(cond_data,2), cond_semi_inds)
        cond_data[:, cond_semi_inds] = convert(eltype(cond_data), NaN)
    end

end

"""
```
string_nearest{T<:AbstactString}(candidates::Vector{T}, target::T)
```
Returns the vector of strings in `candidates` that have the minimum Levenshtein distance
from `target`.
"""
function string_nearest{T<:AbstractString}(candidates::Vector{T}, target::T)
    nearest_str = []
    nearest_dist = typemax(Int)
    for k in keys(candidates)
        this_str = string(k)
        this_dist = levenshtein(this_str, target)
        if this_dist < nearest_dist
            nearest_str = [this_str]
            nearest_dist = this_dist
        elseif this_dist == nearest_dist
            push!(nearest_str, this_str)
        end
    end
    return nearest_str
end
