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
    s = replace(s, "α", "alpha")
    s = replace(s, "β", "beta")
    s = replace(s, "γ", "gamma")
    s = replace(s, "δ", "delta")
    s = replace(s, "ϵ", "epsilon")
    s = replace(s, "ε", "epsilon")
    s = replace(s, "η", "eta")
    s = replace(s, "θ", "theta")
    s = replace(s, "ι", "iota")
    s = replace(s, "κ", "kappa")
    s = replace(s, "λ", "lambda")
    s = replace(s, "μ", "mu")
    s = replace(s, "ν", "nu")
    s = replace(s, "ξ", "xi")
    s = replace(s, "π", "pi")
    s = replace(s, "ρ", "rho")
    s = replace(s, "σ", "sigma")
    s = replace(s, "τ", "tau")
    s = replace(s, "υ", "upsilon")
    s = replace(s, "ϕ", "phi")
    s = replace(s, "φ", "phi")
    s = replace(s, "χ", "chi")
    s = replace(s, "ψ", "psi")
    s = replace(s, "ω", "omega")

    s = replace(s, "Α", "Alpha")
    s = replace(s, "Β", "Beta")
    s = replace(s, "Γ", "Gamma")
    s = replace(s, "Δ", "Delta")
    s = replace(s, "Ε", "Epsilon")
    s = replace(s, "Η", "Eta")
    s = replace(s, "Θ", "Theta")
    s = replace(s, "Ι", "Iota")
    s = replace(s, "Κ", "Kappa")
    s = replace(s, "Λ", "Lambda")
    s = replace(s, "Μ", "Mu")
    s = replace(s, "Ν", "Nu")
    s = replace(s, "Ξ", "Xi")
    s = replace(s, "Π", "Pi")
    s = replace(s, "Ρ", "Rho")
    s = replace(s, "Σ", "Sigma")
    s = replace(s, "Τ", "Tau")
    s = replace(s, "Υ", "Upsilon")
    s = replace(s, "Φ", "Phi")
    s = replace(s, "Χ", "Chi")
    s = replace(s, "Ψ", "Psi")
    s = replace(s, "Ω", "Omega")

    return s
end

function detexify(s::ASCIIString)
    s
end

function detexify(s::Symbol)
    symbol(detexify(string(s)))
end