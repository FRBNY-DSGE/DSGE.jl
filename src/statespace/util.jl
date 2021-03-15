# Some utility functions to keep code clean in statespace_functions.jl

function compute_gensys_gensys2_regimes(m::AbstractDSGEModel)

    regime_switching = haskey(get_settings(m), :regime_switching) && get_setting(m, :regime_switching)
    n_regimes        = (regime_switching && haskey(get_settings(m), :n_regimes)) ? get_setting(m, :n_regimes) : 1

    if haskey(get_settings(m), :gensys2) && get_setting(m, :gensys2)
        if haskey(get_settings(m), :regime_eqcond_info)
            sorted_eqcond = sort!(collect(get_setting(m, :regime_eqcond_info)), by = x -> x[1]) # x -> x[1] sorts by keys
            first_gensys2_ind = findfirst(x -> x[2].alternative_policy.key in # finds the first value (EqcondEntry) whose key
                                          get_setting(m, :temporary_altpolicy_names), sorted_eqcond) # is a temporary altpolicy
            first_gensys2_regime = !isnothing(first_gensys2_ind) ? sorted_eqcond[first_gensys2_ind][1] : nothing
        else
            first_gensys2_regime = nothing
        end

        if isnothing(first_gensys2_regime)
            # throw(GensysError("No equilibrium conditions in any regime are being temporarily replaced, " *
            #                   "but the setting :gensys2 is true."))
            gensys2_regimes = Vector{UnitRange{Int}}(undef, 0)
            gensys_regimes  = UnitRange{Int}[1:n_regimes]
        else
            last_gensys2_regime = haskey(get_settings(m), :temporary_altpolicy_length) ?
                min(first_gensys2_regime + get_setting(m, :temporary_altpolicy_length), n_regimes) :
                n_regimes # NOTE removed a +1 here--if tests start failing, check here first

            gensys_regimes = UnitRange{Int}[1:(first_gensys2_regime - 1)]
            if last_gensys2_regime != n_regimes
                append!(gensys_regimes, [(last_gensys2_regime + 1):n_regimes])
            end
            gensys2_regimes = [first_gensys2_regime-1:last_gensys2_regime]
        end
    else
        gensys2_regimes = Vector{UnitRange{Int}}(undef, 0)
        gensys_regimes  = UnitRange{Int}[1:n_regimes]
    end

    return gensys_regimes, gensys2_regimes
end

"""
```
compute_system_function(system::System{S}) where S<:AbstractFloat
```

### Inputs

- `system::System`: The output of compute_system(m), i.e. the matrix outputs from solving a given model, m.

### Outputs

- `Φ::Function`: transition equation
- `Ψ::Function`: measurement equation
- `F_ϵ::Distributions.MvNormal`: shock distribution
- `F_u::Distributions.MvNormal`: measurement error distribution
"""
function compute_system_function(system::System{S}) where S<:AbstractFloat
    # Unpack system
    TTT    = system[:TTT]
    RRR    = system[:RRR]
    CCC    = system[:CCC]
    QQ     = system[:QQ]
    ZZ     = system[:ZZ]
    DD     = system[:DD]
    EE     = system[:EE]

    # Define transition and measurement functions
    @inline Φ(s_t1::Vector{S}, ϵ_t::Vector{S}) = TTT*s_t1 + RRR*ϵ_t + CCC
    @inline Ψ(s_t::Vector{S}) = ZZ*s_t + DD

    # Define shock and measurement error distributions
    nshocks = size(QQ, 1)
    nobs    = size(EE, 1)
    F_ϵ = Distributions.MvNormal(zeros(nshocks), QQ)
    F_u = Distributions.MvNormal(zeros(nobs),    EE)

    return Φ, Ψ, F_ϵ, F_u
end

function zero_system_constants(system::System{S}) where S<:AbstractFloat
    system = copy(system)

    system.transition.CCC = zeros(size(system[:CCC]))
    system.measurement.DD = zeros(size(system[:DD]))
    system.pseudo_measurement.DD_pseudo = zeros(size(system[:DD_pseudo]))

    return system
end

function zero_system_constants(system::RegimeSwitchingSystem{S}) where S<:AbstractFloat
    system = copy(system)

    for i in 1:n_regimes(system)
        system.transitions[i].CCC = zeros(size(system[i, :CCC]))
        system.measurements[i].DD = zeros(size(system[i, :DD]))
        system.pseudo_measurements[i].DD_pseudo = zeros(size(system[i, :DD_pseudo]))
    end

    return system
end

"""
```
k_periods_ahead_expectations(TTT, CCC, TTTs, CCCs, t, k, permanent_t = length(TTTs))
```

calculates the matrices associated with the expected state `k` periods ahead from `t`.
This function should NOT be used with linear state space system matrices with any unit roots.

The `TTT` and `CCC` inputs are the transition matrix and constant vector associated with
the current period `t`, while the `TTTs` and `CCCs` are vectors containing the time-varying
transition matrices and constant vectors, such that `TTTs[t]` retrieves the time-varying
transition matrix associated with period `t` and `TTTs[t + k]` retrieves the time-varying
transition matrix associated with period `t + k`. The optional argument `permanent_t`
indicates the period for which the matrices/vectors are no longer time-varying, i.e.
if `t >= permanent_t`, then `TTTs[permanent_t]` is the transition matrix.

The formula implemented by this function is
```
𝔼ₜ[sₜ₊ₖ] = (∏ⱼ=₁ᵏ Tₜ₊ⱼ) sₜ + (∑ₘ₌₁ᵏ⁻¹ (∏ⱼ₌ₘ₊₁ᵏ Tₜ₊ⱼ) Cₜ₊ₘ) + Cₜ₊ₖ.
```
Additional simplifications are made if it is known that `t + k > permanent_t`
since this implies some matrices are the same. This recognition reduces
unnecessary computations.
"""
function k_periods_ahead_expectations(TTT::AbstractMatrix, CCC::AbstractVector,
                                      TTTs::Vector{<: AbstractMatrix}, CCCs::Vector{<: AbstractVector},
                                      t::Int, k::Int, permanent_t::Int = length(TTTs);
                                      integ_series::Bool = false)

    if isempty(TTTs) || isempty(CCCs)
        Tᵏ = TTT^k
        if all(CCC .≈ 0.)
            return Tᵏ, CCC
        else
            if integ_series
                T_memo = Dict{Int, typeof(TTT)}()
                T_memo[1] = TTT
                for i in 2:k
                    T_memo[i] = T_memo[i - 1] * TTT
                end
                Tᵏsum = sum([T_memo[j] for j in 1:k])
            else
                Tᵏsum = (I - TTT) \ (I - Tᵏ)
            end
            return Tᵏ, Tᵏsum * CCC
        end
    else
        if t + k <= permanent_t
            # Cannot save computation speed by not calculating further times b/c always time-varying
            T_memo = Dict{Int, eltype(TTTs)}()
            T_memo[k] = TTTs[t + k]
            for i in (k-1):-1:1
                T_memo[i] = T_memo[i + 1] * TTTs[t + i]
            end

            C_accum = deepcopy(CCCs[t + k])
            for i in 1:(k - 1)
                C_accum .+= T_memo[i + 1] * CCCs[t + i]
            end

            return T_memo[1], C_accum
        elseif t == permanent_t
            # None of the matrices are time-varying anymore
            if all(CCCs[permanent_t] .≈ 0.)
                return TTTs[permanent_t]^k, CCCs[permanent_t]
            else
                Tᵏₜ₊₁ = TTTs[permanent_t]^k
                if integ_series
                    T_memo = Dict{Int, eltype(TTTs)}()
                    T_memo[1] = TTTs[permanent_t]
                    for i in 2:k
                        T_memo[i] = T_memo[i - 1] * TTTs[permanent_t]
                    end
                    Tᵏsum = sum([T_memo[j] for j in 1:k])
                else
                    Tᵏsum = (I - TTTs[permanent_t]) \ (I - Tᵏₜ₊₁)
                end
                return Tᵏₜ₊₁, Tᵏsum * CCCs[permanent_t]
            end
        else
            # Computation time can be saved by realizing some matrices are not time-varying
            h = (permanent_t - 1) - t # last time of time-variation is permanent_t - 1
            Tᵏ⁻ʰₜ₊ₕ₊₁ = (k == h) ? Diagonal(ones(length(CCCs[permanent_t]))) : TTTs[permanent_t]^(k - h)

            T_memo = Dict{Int, eltype(TTTs)}()
            if h > 0
                T_memo[h] = TTTs[t + h] # maps i to ∏ⱼ₌ᵢʰ Tₜ₊ⱼ, so T_memo[h] = Tₜ₊ₕ, T_memo[h-1] = Tₜ₊ₕ * Tₜ₊ₕ₋₁, ...
                for i in (h-1):-1:1
                    T_memo[i] = T_memo[i + 1] * TTTs[t + i]
                end
                T_accum = Tᵏ⁻ʰₜ₊ₕ₊₁ * T_memo[1]
            else
                T_accum = Tᵏ⁻ʰₜ₊ₕ₊₁
            end

            if h == 0 # Nothing to accumulate from the past
                C_accum = zeros(size(T_accum, 1))
            else
                C_accum = deepcopy(CCCs[t + h])
                for i in 1:(h - 1)
                    C_accum .+= T_memo[i + 1] * CCCs[t + i]
                end
                C_accum .= Tᵏ⁻ʰₜ₊ₕ₊₁ * C_accum
            end

            if all(CCCs[permanent_t] .≈ 0.)
                return T_accum, C_accum
            else
                if integ_series
                    T_memo = Dict{Int, eltype(TTTs)}()
                    T_memo[1] = TTTs[permanent_t]
                    for i in 2:(k - h - 1)
                        T_memo[i] = T_memo[i - 1] * TTTs[permanent_t]
                    end
                    Tᵏ⁻ᵐsum = sum([T_memo[j] for j in 1:(k - h - 1)]) + I
                else
                    Tᵏ⁻ᵐsum = (I - TTTs[permanent_t]) \ (I - Tᵏ⁻ʰₜ₊ₕ₊₁)
                end
                C_accum .+= Tᵏ⁻ᵐsum * CCCs[permanent_t]
                return T_accum, C_accum
            end
        end
    end
end

"""
```
k_periods_ahead_expected_sums(TTT, CCC, TTTs, CCCs, t, k, permanent_t = length(TTTs))
```

calculates the matrices associated with the sum of the expected states between periods
`t + 1` and  `t + k`. This function should NOT be used with
linear state space system matrices with any unit roots.

The `TTT` and `CCC` inputs are the transition matrix and constant vector associated with
the current period `t`, while the `TTTs` and `CCCs` are vectors containing the time-varying
transition matrices and constant vectors, such that `TTTs[t]` retrieves the time-varying
transition matrix associated with period `t` and `TTTs[t + k]` retrieves the time-varying
transition matrix associated with period `t + k`. The optional argument `permanent_t`
indicates the period for which the matrices/vectors are no longer time-varying, i.e.
if `t >= permanent_t`, then `TTTs[permanent_t]` is the transition matrix.

The formula implemented by this function is
```
∑ⱼ₌₁ᵏ 𝔼ₜ[sₜ₊ⱼ] = ∑ⱼ₌₁ᵏ(∏ⱼ=₁ᵏ Tₜ₊ⱼ) sₜ + ∑ᵣ₌₁ᵏ⁻¹(I + ∑ⱼ₌ᵣ₊₁ᵏ (∏ₘ₌ᵣ₊₁ʲ Tₜ₊ₘ))Cₜ₊ᵣ + Cₜ₊ₖ.
```
Additional simplifications are made if it is known that `t + k > permanent_t`
since this implies some matrices are the same. This recognition reduces
unnecessary computations.
"""
function k_periods_ahead_expected_sums(TTT::AbstractMatrix, CCC::AbstractVector,
                                       TTTs::Vector{<: AbstractMatrix}, CCCs::Vector{<: AbstractVector},
                                       t::Int, k::Int, permanent_t::Int = length(TTTs);
                                       integ_series::Bool = false)

    if integ_series # Do this by directly summing the k-periods ahead expectations. Not fully efficient but also not the typical case
        T_accum = Vector{eltype(TTTs)}(undef, k)
        C_accum = Vector{eltype(CCCs)}(undef, k)
        for i in 1:k
            T_accum[i], C_accum[i] = k_periods_ahead_expectations(TTT, CCC, TTTs, CCCs, t, i,
                                                                  integ_series = integ_series)
        end
        return sum(T_accum), sum(C_accum)
    end

    if isempty(TTTs) || isempty(CCCs)
        Tᵏsum = (I - TTT) \ (TTT - TTT^(k + 1))
        if all(CCC .≈ 0.)
            return Tᵏsum, CCC
        else
            TTTʲmemo = Vector{typeof(TTT)}(undef, k)
            TTTʲmemo[1] = TTT
            for q in 2:k
                TTTʲmemo[q] = TTTʲmemo[q - 1] * TTT
            end
            return Tᵏsum, (I + sum([(I - TTT) \ (I - TTTʲmemo[k - q + 1]) for q in 1:(k - 1)])) * CCC
        end
    else
        if t + k <= permanent_t
            # Cannot save computation speed by not calculating further times b/c always time-varying
            total_Tsum = zeros(eltype(TTTs[t]), size(TTTs[t]))
            total_Csum = zeros(eltype(CCCs[t]), size(CCCs[t]))

            T_accum = copy(TTTs[t+k])
            C_accum = copy(CCCs[t+k])
            total_Tsum .+= T_accum
            total_Csum .+= C_accum

            for i in (k-1):-1:1
                C_accum .+= T_accum * CCCs[t+i]
                T_accum .*= TTTs[t+i]

                total_Tsum .+= T_accum
                total_Csum .+= C_accum
            end

            return total_Tsum, total_Csum
        elseif t == permanent_t
            # None of the matrices are time-varying anymore
            Tᵏsum = (I - TTTs[permanent_t]) \ (TTTs[permanent_t] - TTTs[permanent_t]^(k + 1))
            if all(CCC .≈ 0.)
                return Tᵏsum, CCCs[permanent_t]
            else
                TTTʲmemo = Vector{eltype(TTTs)}(undef, k)
                TTTʲmemo[1] = TTTs[permanent_t]
                for q in 2:k
                    TTTʲmemo[q] = TTTʲmemo[q - 1] * TTTs[permanent_t]
                end
                return Tᵏsum, (I + sum([(I - TTTs[permanent_t]) \ (I - TTTʲmemo[k - q + 1]) for q in 1:(k - 1)])) * CCCs[permanent_t]
            end
        else
            total_Tsum = zeros(eltype(TTTs[t]), size(TTTs[t]))
            total_Csum = zeros(eltype(CCCs[t]), size(CCCs[t]))

            # T_memo and C_memo store the expectation of T and C at time i
            ## T_memo[i] will store ∏_{j=t+1}^i T_j
            ## C_memo[i] stores ∑_{q=t+1}^i C_q + ∑_{j=t+1}^q T_j C_q
            # total_sum will stores sum over the expectations
            T_memo = Dict{Int, eltype(TTTs)}()
            C_memo = Dict{Int, eltype(CCCs)}()
            T_memo[t+1] = copy(TTTs[t+1])
            C_memo[t+1] = copy(CCCs[t+1])
            total_Tsum .+= T_memo[t+1]
            total_Csum .+= C_memo[t+1]

            do_C = !(all(CCCs[permanent_t] .≈ 0.0))

            for i in (t+2):permanent_t
                T_memo[i] = TTTs[i] * T_memo[i-1]
                total_Tsum .+= T_memo[i]

                C_memo[i] = TTTs[i] * C_memo[i-1] + CCCs[i]
                total_Csum .+= C_memo[i]
            end

            T_sums = (I - TTTs[permanent_t]) \ (TTTs[permanent_t] - TTTs[permanent_t]^(t+k-permanent_t+1))
            total_Tsum .+= T_sums * T_memo[permanent_t]

            if !do_C
                C_tmp = T_sums * C_memo[permanent_t]
                total_Csum .+= C_tmp
            else
                for i in (permanent_t+1):(t+k)
                    C_memo[i] = TTTs[permanent_t] * C_memo[i-1] + CCCs[permanent_t]
                    total_Csum .+= C_memo[i]
                end
            end
#=
            if t == 6
                # Testing code
                T_accum2 = Vector{eltype(TTTs)}(undef, k)
                C_accum2 = Vector{eltype(CCCs)}(undef, k)
                for i in 1:k
                    T_accum2[i], C_accum2[i] = k_periods_ahead_expectations(TTT, CCC, TTTs, CCCs, t, i,
                                                                          integ_series = integ_series)
                end
                JLD2.jldopen("testing_accum_full.jld2", true, true, true, IOStream) do file
                    write(file, "T_accum", T_accum2)
                    write(file, "C_accum", C_accum2)
                    write(file, "T_new", T_memo)
                    write(file, "C_new", C_memo)
                    write(file, "TTTs", TTTs)
                    write(file, "CCCs", CCCs)
                    write(file, "T_sum", total_Tsum)
                    write(file, "C_sum", total_Csum)
                end
            end
=#
            return total_Tsum, total_Csum

            # This code seems to work (match directly adding each k_period_ahead_expectations)
            # but sometimes causes mild numerical differences,
            # so we calculate the expected sum by adding each k_period_ahead_expectations.
            #=            h = (permanent_t - 1) - t # last time of time-variation is permanent_t - 1

            Tₜ₊ₕ₊₁_memo = Vector{eltype(TTTs)}(undef, k - h)
            Tₜ₊ₕ₊₁_memo[1] = TTTs[permanent_t] # maps j to Tₜ₊ₕ₊₁ʲ⁻ʰ for j in (h + 1):k or Tₜ₊ₕ₊₁ʲ for j in 1:(k-h)
            for j in 2:(k - h)
            Tₜ₊ₕ₊₁_memo[j] = TTTs[permanent_t] * Tₜ₊ₕ₊₁_memo[j - 1]
            end

            Tₜ₊ₕ₊₁ʲsum = (I - TTTs[permanent_t]) \ (TTTs[permanent_t] - TTTs[permanent_t] ^ (k - h + 1)) # ∑ⱼ₌₁ᵏ⁻ʰ Tₜ₊ₕ₊₁ʲ

            TC_memo = Vector{eltype(TTTs)}(undef, h) # matrices to be used for both accumulated TTT and CCC
            TC_memo[h] = TTTs[t + h] # maps i to ∏ⱼ₌ᵢʰ Tₜ₊ⱼ, so T_memo[h] = Tₜ₊ₕ, T_memo[h-1] = Tₜ₊ₕ * Tₜ₊ₕ₋₁, ...
            for i in (h-1):-1:1
            TC_memo[i] = TC_memo[i + 1] * TTTs[t + i]
            end

            T1_term = Vector{eltype(TTTs)}(undef, h) # first part of the accumulated TTT matrix
            T1_term[1] = TTTs[t + 1] # maps i to ∏ⱼ₌₁ⁱ Tₜ₊ⱼ, so T_memo[h] = Tₜ₊ₕ * ⋯  * Tₜ₊₁
            for i in 2:(h - 1)
            T1_term[i] = TTTs[t + i] * T1_term[i - 1]
            end
            T1_term[h] = TC_memo[1] # This one was calculated already

            # second part of the accumulated TTT matrix
            # maps j to Tₜ₊ₕ₊₁ʲ ∏ᵣ₌₁ʰ Tₜ₊ᵣ, so T_memo[h] = Tₜ₊ₕ * ⋯  * Tₜ₊₁
            T2_term = [Tₜ₊ₕ₊₁_memo[j] * T1_term[h] for j in 1:(k - h)]

            T_accum = sum(T1_term) + sum(T2_term) # calculated accumulated matrix

            # Calculate final portions of accumulated CCC vector first
            I₊Tʲ = (I + Tₜ₊ₕ₊₁ʲsum)
            C_accum   = I₊Tʲ * CCCs[t + h]
            if any(.!(CCCs[permanent_t] .≈ 0.))
            C_accum .+= (I + (k - h - 1) .* I₊Tʲ) * CCCs[permanent_t]
            end
            if h > 1
            C1_term = Vector{eltype(TTTs)}(undef, h - 1)
            C2_term = Vector{eltype(TTTs)}(undef, h - 1)

            for q in 1:(h - 1)
            C1_term[q] = sum([prod([TTTs[t + m] for m in (q + 1):j]) for j in (q + 1):h])
            C2_term[q] = sum([Tₜ₊ₕ₊₁_memo[j] * TC_memo[q + 1] for j in 1:(k - h)])
            end
            C_accum .+= sum([(I + C1_term[q] + C2_term[q]) * CCCs[t + q] for q in 1:(h - 1)])
            end=#
        end
    end
end
