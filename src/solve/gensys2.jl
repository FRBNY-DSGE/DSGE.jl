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

    for row in 1:length(exp_eq_ind)
        if exp_eq_ind[row] == 0
            # Not expectational equation
            Γ1_til[row, :] = Γ1[row, :]
            Γ0_til[row, :] = Γ0[row, :]
            Γ2_til[row, :] .= 0.
        else
            # expectational equation
            Γ0_til[row, findfirst(Γ1[row, :] .> 0)] = -1
            Γ2_til[row, findfirst(Γ0[row, :] .> 0)] = 1
        end
    end

#     Γ0_til, Γ1_til, Γ2_til, C_til, Ψ_til = gensys_to_predictable_form(Γ0, Γ1, C, Ψ, Π)

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

    Γ0_tils = Vector{Matrix{Float64}}(undef, length(Γ0s))
    Γ1_tils = Vector{Matrix{Float64}}(undef, length(Γ0s))
    Γ2_tils = Vector{Matrix{Float64}}(undef, length(Γ0s))
    C_tils = Vector{Vector{Float64}}(undef, length(Γ0s))
    Ψ_tils = Vector{Matrix{Float64}}(undef, length(Γ0s))

    for i in 1:length(Γ0s)
        exp_eq_ind = sum(Πs[i], dims = 2)
        Γ0_tils[i] = zeros(size(Γ0s[i]))
        Γ1_tils[i] = zeros(size(Γ1s[i]))
        Γ2_tils[i] = zeros(size(Γ0s[i]))
        C_tils[i] = Cs[i]
        Ψ_tils[i] = Ψs[i]

        for row in 1:length(exp_eq_ind)
            if exp_eq_ind[row] == 0
                # Not expectational equation
                Γ1_tils[i][row, :] = Γ1s[i][row, :]
                Γ0_tils[i][row, :] = Γ0s[i][row, :]
                Γ2_tils[i][row, :] .= 0.
            else
                # expectational equation
                Γ0_tils[i][row, findfirst(Γ1s[i][row, :] .> 0)] = -1
                Γ2_tils[i][row, findfirst(Γ0s[i][row, :] .> 0)] = 1
            end
        end
    end

#=    for i in 1:length(Γ0s)
        Γ0_tils[i], Γ1_tils[i], Γ2_tils[i], C_tils[i], Ψ_tils[i] = gensys_to_predictable_form(Γ0s[i], Γ1s[i], Cs[i], Ψs[i], Πs[i])
    end=#

    Tcal = Vector{Matrix{Float64}}(undef, T_switch)
    Rcal = Vector{Matrix{Float64}}(undef, T_switch)
    Ccal = Vector{Vector{Float64}}(undef, T_switch)

    Tcal[end] = TTT
    Rcal[end] = RRR
    Ccal[end] = CCC

    for t = 1:(T_switch-1)
        Tcal[end-t] = (Γ2_tils[end-t]*TTT + Γ0_tils[end-t])\Γ1_tils[end-t]
        Rcal[end-t] = (Γ2_tils[end-t]*TTT + Γ0_tils[end-t])\Ψ_tils[end-t]
        Ccal[end-t] = (Γ2_tils[end-t]*TTT + Γ0_tils[end-t])\(C_tils[end-t] - Γ2_tils[end-t]*CCC)

        TTT = Tcal[end-t]
        CCC = Ccal[end-t]
        RRR = Rcal[end-t]
    end

    return Tcal, Rcal, Ccal
end

function gensys_to_predictable_form(Γ0::AbstractMatrix{S}, Γ1::AbstractMatrix{S}, C::AbstractVector{S},
                                    Ψ::AbstractMatrix{S}, Π::AbstractMatrix{S}) where {S <: Real}

    nonexp_eq_ind = map(i -> all(Π[i, :] .== 0.), 1:size(Π, 1))
    Γ0_til = zeros(S, size(Γ0))
    Γ1_til = zeros(S, size(Γ1))
    Γ2_til = zeros(S, size(Γ0))
    C_til  = C
    Ψ_til  = Ψ

    for row in 1:length(nonexp_eq_ind)
        if nonexp_eq_ind[row]
            # Not expectational equation
            Γ1_til[row, :]  = Γ1[row, :]
            Γ0_til[row, :]  = Γ0[row, :]
        else
            # expectational equation, assumed to take the form zₜ = Eₜ₋₁[zₜ] + ηₜ
            Γ0_til[row, findfirst(Γ1[row, :] .> 0)] = -1.
            Γ2_til[row, findfirst(Γ0[row, :] .> 0)] = 1.
        end
    end

    return Γ0_til, Γ1_til, Γ2_til, C_til, Ψ_til
end
