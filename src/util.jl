"""
```
sendto(p::Int; args...)
```
Function to send data from master process to particular worker, p. Code from ChrisRackauckas, avavailable at: https://github.com/ChrisRackauckas/ParallelDataTransfer.jl/blob/master/src/ParallelDataTransfer.jl.
"""
function sendto(p::Int; args...)
    for (nm, val) in args
        @spawnat(p, Core.eval(Main, Expr(:(=), nm, val)))
    end
end

"""
```
sendto(ps::AbstractVector{Int}; args...)
```
Function to send data from master process to list of workers. Code from ChrisRackauckas, available at: https://github.com/ChrisRackauckas/ParallelDataTransfer.jl/blob/master/src/ParallelDataTransfer.jl.
"""
function sendto(ps::AbstractVector{Int}; args...)
    for p in ps
        sendto(p; args...)
    end
end

"""
```
abbrev_symbol(s::Symbol, n::Int=4)
```

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

function sorted_list_insert!(v::Vector{T}, x::T) where T
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
    dr = t0:Day(1):t1
    return Base.filter(d -> Dates.lastdayofquarter(d) == d, dr)
end

"""
```
detexify(s::String)

detexify(s::Symbol)
```

Remove Unicode characters from the string `s`, replacing them with ASCII
equivalents. For example, `detexify(\"π\")` returns `\"pi\"`.
"""
function detexify(s::String)
    s = replace(s, "α" => "alpha")
    s = replace(s, "β" => "beta")
    s = replace(s, "γ" => "gamma")
    s = replace(s, "δ" => "delta")
    s = replace(s, "ϵ" => "epsilon")
    s = replace(s, "ε" => "epsilon")
    s = replace(s, "ζ" => "zeta")
    s = replace(s, "η" => "eta")
    s = replace(s, "θ" => "theta")
    s = replace(s, "ι" => "iota")
    s = replace(s, "κ" => "kappa")
    s = replace(s, "λ" => "lambda")
    s = replace(s, "μ" => "mu")
    s = replace(s, "ν" => "nu")
    s = replace(s, "ξ" => "xi")
    s = replace(s, "π" => "pi")
    s = replace(s, "ρ" => "rho")
    s = replace(s, "ϱ" => "rho")
    s = replace(s, "σ" => "sigma")
    s = replace(s, "ς" => "sigma")
    s = replace(s, "τ" => "tau")
    s = replace(s, "υ" => "upsilon")
    s = replace(s, "ϕ" => "phi")
    s = replace(s, "φ" => "varphi")
    s = replace(s, "χ" => "chi")
    s = replace(s, "ψ" => "psi")
    s = replace(s, "ω" => "omega")

    s = replace(s, "Α" => "Alpha")
    s = replace(s, "Β" => "Beta")
    s = replace(s, "Γ" => "Gamma")
    s = replace(s, "Δ" => "Delta")
    s = replace(s, "Ε" => "Epsilon")
    s = replace(s, "Ζ" => "Zeta")
    s = replace(s, "Η" => "Eta")
    s = replace(s, "Θ" => "Theta")
    s = replace(s, "Ι" => "Iota")
    s = replace(s, "Κ" => "Kappa")
    s = replace(s, "Λ" => "Lambda")
    s = replace(s, "Μ" => "Mu")
    s = replace(s, "Ν" => "Nu")
    s = replace(s, "Ξ" => "Xi")
    s = replace(s, "Π" => "Pi")
    s = replace(s, "Ρ" => "Rho")
    s = replace(s, "Σ" => "Sigma")
    s = replace(s, "Τ" => "Tau")
    s = replace(s, "Υ" => "Upsilon")
    s = replace(s, "Φ" => "Phi")
    s = replace(s, "Χ" => "Chi")
    s = replace(s, "Ψ" => "Psi")
    s = replace(s, "Ω" => "Omega")

    s = replace(s, "′" => "'")

    return s
end

function detexify(s::Symbol)
    Symbol(detexify(string(s)))
end

dispfns = [:print, :println]
for disp in dispfns
    @eval begin
        function Base.$disp(verbose::Symbol, min::Symbol, xs...)
            if VERBOSITY[verbose] >= VERBOSITY[min]
                $disp(xs...)
            end
        end
    end
end

function info_print(verbose::Symbol, min::Symbol, x::String)
    if VERBOSITY[verbose] >= VERBOSITY[min]
        @info(x)
    end
end

function warn_print(verbose::Symbol, min::Symbol, x::String)
    if VERBOSITY[verbose] >= VERBOSITY[min]
        @warn(x)
    end
end

## Testing functions

function test_matrix_eq2(expect::AbstractArray,
                         actual::AbstractArray,
                         expstr::String,
                         actstr::String,
                         ϵ_abs::Float64 = 1e-6,
                         ϵ_rel::Float64 = 1e-2) where {T<:AbstractFloat}
    if length(expect) ≠ length(actual)
        error("") #="lengths of ", expstr, " and ", actstr, " do not match: ",
              "\n  ", expstr, " (length $(length(expect))) = ", expect,
              "\n  ", actstr, " (length $(length(actual))) = ", actual)=#
    end

    # Absolute difference filter
    abs_diff   = abs.(actual .- expect) .> ϵ_abs
    n_abs_diff = sum(skipmissing(abs_diff))

    # Relative difference filter
    rel_diff   = 100*abs.((actual .- expect) ./ expect) .> ϵ_rel
    n_rel_diff = sum(skipmissing(rel_diff))

    # Element is only problematic if it fails *both* tests.
    mixed_diff   = abs_diff .& rel_diff
    n_mixed_diff = sum(skipmissing(mixed_diff))

    if n_mixed_diff ≠ 0
        sdiff = string("|a - b| <= ", ϵ_abs,
                   " or |a - b|/|b| <= ", ϵ_rel, "%,",
                   " ∀ a ∈ ", actstr, ",",
                   " ∀ b ∈ ",expstr)
        @warn "assertion failed:\n",
             "    ", sdiff,
             "\n$(n_abs_diff) entries fail absolute equality filter",
             "\n$(n_rel_diff) entries fail relative equality filter",
             "\n$(n_mixed_diff) entries fail both equality filters\n"
        return false
    end
    return true
end

"""
    @test_matrix_approx_eq(a, b)

Test two matrices of floating point numbers `a` and `b` for approximate equality.
"""
macro test_matrix_approx_eq(a,b)
    :(test_matrix_eq2($(esc(a)),$(esc(b)),$(string(a)),$(string(b))))
end

"""
    @test_matrix_approx_eq_eps(a, b, ϵ_abs, ϵ_rel)

Test two matrices of floating point numbers `a` and `b` for equality taking in account a
margin of absolute tolerance given by `ϵ_abs` and a margin of relative tolerance given by
`ϵ_rel` (in comparison to `b`).
"""
macro test_matrix_approx_eq_eps(a,b,c,d)
    :(test_matrix_eq2($(esc(a)),$(esc(b)),$(string(a)),$(string(b)),$(esc(c)),$(esc(d))))
end

"""
Sparse identity matrix - since deprecated in 0.7
"""
function speye(n::Integer)
    return SparseMatrixCSC{Float64}(I, n, n)
end

"""
Sparse identity matrix - since deprecated in 0.7
"""
function speye(T::Type, n::Integer)
    return SparseMatrixCSC{T}(I, n, n)
end


"""
    <(a::Complex, b::Complex)

Compare real values of complex numbers.
"""
function <(a::Complex, b::Complex)
    return a.re < b.re
end

"""
    <(a::Real, b::Complex)

Compare real values of complex numbers.
"""
function <(a::Real, b::Complex)
    return a < b.re
end

"""
    <(a::Complex, b::Real)

Compare real values of complex numbers.
"""
function <(a::Complex, b::Real)
    return a.re < b
end

function min(a::Complex, b::Real)
    return min(a.re, b)
end

function min(a::Complex, b::Complex)
    return min(a.re, b.re)
end

function min(a::Real, b::Complex)
    return min(a, b.re)
end

function max(a::Complex, b::Real)
    return max(a.re, b)
end

function max(a::Complex, b::Complex)
    return max(a.re, b.re)
end

function max(a::Real, b::Complex)
    return max(a, b.re)
end

"""
```
n_param_regs(params::ParameterVector)
```
Get total number of parameter regimes for each parameter
"""
function n_param_regs(params::ParameterVector)
    return [haskey(params[i].regimes, :value) ? length(params[i].regimes[:value]) : 1 for i in 1:length(m.parameters)]
end
