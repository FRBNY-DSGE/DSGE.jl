"""
```
gensys2(m::AbstractDSGEModel, Γ0::Matrix{Float64}, Γ1::Matrix{Float64}, C::Vector{Float64},
        Ψ::Matrix{Float64}, Π::Matrix{Float64}, TTT::Matrix{Float64}, RRR::Matrix{Float64},
        CCC::Vector{Float64}, T_switch::Int)
```

calculates the state space transition matrices when a temporary alternative policy applies.
This function is **not** the same as the gensys2 written by Chris Sims, which computes
a second-order perturbation.
"""
function gensys2(m::AbstractDSGEModel, Γ0::Matrix{Float64}, Γ1::Matrix{Float64}, C::Vector{Float64},
                 Ψ::Matrix{Float64}, Π::Matrix{Float64}, TTT::Matrix{Float64}, RRR::Matrix{Float64},
                 CCC::Vector{Float64}, T_switch::Int)

    exp_eq_ind = sum(Π, dims = 2)
    Γ0_til = zeros(size(Γ0))
    Γ1_til = zeros(size(Γ1))
    Γ2_til = zeros(size(Γ0))
    C_til = C
    Ψ_til = Ψ

    Γ0_til, Γ1_til, Γ2_til, C_til, Ψ_til =
        gensys_to_predictable_form(Γ0, Γ1, C, Ψ, Π; use_sparse = (haskey(get_settings(m), :gensys2_sparse_matrices) &&
                                                                  get_setting(m, :gensys2_sparse_matrices)))

    Tcal = Vector{Matrix{Float64}}(undef, T_switch)
    Rcal = Vector{Matrix{Float64}}(undef, T_switch)
    Ccal = Vector{Vector{Float64}}(undef, T_switch)

    Tcal[end] = TTT
    Rcal[end] = RRR
    Ccal[end] = CCC

    for t = 1:(T_switch-1)
        Tcal[end-t] = (Γ2_til*TTT + Γ0_til)\Γ1_til
        Rcal[end-t] = (Γ2_til*TTT + Γ0_til)\Ψ_til
        Ccal[end-t] = (Γ2_til*TTT + Γ0_til)\(C_til - Γ2_til*CCC)

        TTT = Tcal[end-t]
        CCC = Ccal[end-t]
        RRR = Rcal[end-t]
    end
    return Tcal, Rcal, Ccal
end

function gensys2(m::AbstractDSGEModel, Γ0s::Vector{Matrix{Float64}}, Γ1s::Vector{Matrix{Float64}},
                 Cs::Vector{Vector{Float64}}, Ψs::Vector{Matrix{Float64}}, Πs::Vector{Matrix{Float64}},
                 TTT::Matrix{Float64}, RRR::Matrix{Float64},
                 CCC::Vector{Float64}, T_switch::Int)

    ntil = length(Γ0s)
    Γ0_tils = Vector{Matrix{Float64}}(undef, ntil)
    Γ1_tils = Vector{Matrix{Float64}}(undef, ntil)
    Γ2_tils = Vector{Matrix{Float64}}(undef, ntil)
    C_tils = Vector{Vector{Float64}}(undef, ntil)
    Ψ_tils = Vector{Matrix{Float64}}(undef, ntil)

    for i in 1:length(Γ0s)
        Γ0_tils[i], Γ1_tils[i], Γ2_tils[i], C_tils[i], Ψ_tils[i] =
            gensys_to_predictable_form(Γ0s[i], Γ1s[i], Cs[i], Ψs[i], Πs[i];
                                       use_sparse = (haskey(get_settings(m), :gensys2_sparse_matrices) &&
                                                     get_setting(m, :gensys2_sparse_matrices)))
    end

    Tcal = Vector{Matrix{Float64}}(undef, T_switch)
    Rcal = Vector{Matrix{Float64}}(undef, T_switch)
    Ccal = Vector{Vector{Float64}}(undef, T_switch)

    # maybe allocate sparse version of TTT b/c going to do a lot of calculations with it
    # and also make sure to check for essentially zero terms.
    # leave Gamma1_til as a full array. Only make Gamma2, Gamma0 sparse.
    Tcal[end] = TTT
    Rcal[end] = RRR
    Ccal[end] = CCC

    for t = 1:(T_switch-1)
        Tcal[end-t] = (Γ2_tils[end-t] * TTT + Γ0_tils[end-t]) \ Γ1_tils[end-t]
        Rcal[end-t] = (Γ2_tils[end-t] * TTT + Γ0_tils[end-t]) \ Ψ_tils[end-t]
        Ccal[end-t] = (Γ2_tils[end-t] * TTT + Γ0_tils[end-t]) \ (C_tils[end-t] - Γ2_tils[end-t] * CCC)

        TTT = Tcal[end-t]
        CCC = Ccal[end-t]
        RRR = Rcal[end-t]
    end

    return Tcal, Rcal, Ccal
end

function gensys_to_predictable_form(Γ0::AbstractMatrix{S}, Γ1::AbstractMatrix{S}, C::AbstractVector{S},
                                    Ψ::AbstractMatrix{S}, Π::AbstractMatrix{S};
                                    use_sparse::Bool = false) where {S <: Real}

    nonexp_eq_ind = map(i -> all(Π[i, :] .== 0.), 1:size(Π, 1))
    if use_sparse
        Γ0_til = spzeros(S, size(Γ0))
        Γ1_til = zeros(S, size(Γ1))   # this one has to be a dense matrix for left-divide
        Γ2_til = spzeros(S, size(Γ0))
    else
        Γ0_til = zeros(S, size(Γ0))
        Γ1_til = zeros(S, size(Γ1))
        Γ2_til = zeros(S, size(Γ0))
    end
    C_til  = C
    Ψ_til  = Ψ

    for row in 1:length(nonexp_eq_ind)
        if nonexp_eq_ind[row]
            # Not expectational equation
            Γ1_til[row, :]  = view(Γ1, row, :)
            Γ0_til[row, :]  = view(Γ0, row, :)
        else
            # expectational equation, assumed to take the form zₜ = Eₜ₋₁[zₜ] + ηₜ
            Γ0_til[row, findfirst(x -> x > 0., view(Γ1, row, :))] = -1.
            Γ2_til[row, findfirst(x -> x > 0., view(Γ0, row, :))] = 1.
        end
    end

    return Γ0_til, Γ1_til, Γ2_til, C_til, Ψ_til
end
