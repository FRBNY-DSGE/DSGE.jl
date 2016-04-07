"""
    abbrev_symbol(s::Symbol, n::Int=4)

Abbreviate the symbol `s` to an `n`-character string.
"""
function abbrev_symbol(s::Symbol, n::Int=4)
    str = string(s)
    if length(str) ≤ n
        return str
    else
        return str[1:n]
    end
end

"""
```
string_nearest{T<:AbstactString}(candidates::Vector{T}, target::T)
```
Return the vector of strings in `candidates` that have the minimum Levenshtein distance
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

function sorted_list_insert!{T}(v::Vector{T}, x::T)
    insert_index = 1
    for val in v
        if x<val
            break
        else
            insert_index += 1
        end
    end
    insert!(v,insert_index,x)
end

