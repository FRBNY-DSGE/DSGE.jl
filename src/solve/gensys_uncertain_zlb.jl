function gensys_uncertain_zlb(prob_vec::AbstractVector{S}, Th::AbstractMatrix{S}, Ch::AbstractVector{S},
                              Tzlbs::Vector{Matrix{S}}, Rzlbs::Vector{Matrix{S}}, Czlbs::Vector{Vector{S}},
                              Γ0_til::AbstractMatrix{S}, Γ1_til::AbstractMatrix{S}, Γ2_til::AbstractMatrix{S},
                              C_til::AbstractVector{S}, Ψ_til::AbstractMatrix{S}) where {S <: Real}

    Tbars = Vector{Matrix{S}}(undef, length(Tzlbs) + 1)
    Cbars = Vector{Vector{S}}(undef, length(Tzlbs) + 1)
    Tout  = Vector{Matrix{S}}(undef, length(Tzlbs) + 1)
    Rout  = Vector{Matrix{S}}(undef, length(Tzlbs) + 1)
    Cout  = Vector{Vector{S}}(undef, length(Tzlbs) + 1)

    # Calculate "uncertain" ZLB matrices and back out the transition equation
    for i in 1:length(Tzlbs)
        Tbars[i] = prob_vec[1] * Tzlbs[i] + prob_vec[2] * Th # it is assumed Tzlbs is a vector of the T_{t + 1}^{(zlb)} matrices
        Cbars[i] = prob_vec[1] * Czlbs[i] + prob_vec[2] * Ch

        Tout[i]  = (Γ2_til * Tbars[i] + Γ0_til) \ Γ1_til
        Rout[i]  = (Γ2_til * Tbars[i] + Γ0_til) \ Ψ_til
        Cout[i]  = (Γ2_til * Tbars[i] + Γ0_til) \ (C_til  - Γ2_til * Cbars[i])
    end

    # Add boundary condition
    Tout[end] = Tzlbs[end]
    Rout[end] = Rzlbs[end]
    Cout[end] = Czlbs[end]

    return Tout, Rout, Cout
end
