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

"""
```
quarter_range(t0::Date, t1::Date)
```

Returns a vector of `Dates`, consisting of the last days of each quarter between
`t0` and `t1`, inclusive.
"""
function quarter_range(t0::Date, t1::Date)
    dr = t0:t1
    return Dates.recur(d -> Dates.lastdayofquarter(d) == d, dr)
end

function detexify(s::UTF8String)
    replace(s, "π", "pi")
end

function detexify(s::ASCIIString)
    s
end