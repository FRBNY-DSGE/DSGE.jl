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
